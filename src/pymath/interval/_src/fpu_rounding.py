import ctypes
import math
import sys
from decimal import Decimal
from typing import SupportsFloat

if sys.version_info < (3, 9):
    from typing import Callable, Iterator, List, Tuple
else:
    from builtins import list as List, tuple as Tuple
    from collections.abc import Callable, Iterator

SUBNORMAL = math.ldexp(1.0, -1022)
NEAR_ZERO = math.ldexp(1.0, -1074)

nextafter: Callable[[float, float], float]

try:
    if sys.version_info >= (3, 9):
        from math import nextafter
    elif sys.platform == "linux" or sys.platform == "linux2":
        nextafter = ctypes.cdll.LoadLibrary("libm.so.6").nextafter
        nextafter.restype = ctypes.c_double
        nextafter.argtypes = [ctypes.c_double, ctype.c_double]
    elif sys.platform == "darwin":
        nextafter = ctypes.cdll.LoadLibrary("libSystem.dylib").nextafter
        nextafter.restype = ctypes.c_double
        nextafter.argtypes = [ctypes.c_double, ctype.c_double]
    elif sys.platform == "win32":
        nextafter = ctypes.cdll.LoadLibrary("libSystem.dylib")._nextafter
        nextafter.restype = ctypes.c_double
        nextafter.argtypes = [ctypes.c_double, ctype.c_double]
    else:
        from math import nextafter
except:
    def nextafter(x: float, y: float) -> float:
        if isinstance(x, SupportsFloat):
            x = float(x)
        if isinstance(y, SupportsFloat):
            y = float(y)
        if x == y:
            return x
        elif math.isnan(x) or math.isnan(y):
            return x + y
        elif math.isinf(x):
            return x
        elif -SUBNORMAL < x < SUBNORMAL:
            if x < y:
                return x - NEAR_ZERO
            else:
                return x + NEAR_ZERO
        mantissa, exponent = math.frexp(x)
        if x < y:
            mantissa += 0.5 * sys.float_info.epsilon
        else:
            mantissa -= 0.5 * sys.float_info.epsilon
        return math.ldexp(mantissa, exponent)

from .typing import SupportsRichFloat

def split_bits(n: int) -> Iterator[int]:
    if n > 0:
        for i in range(0, n.bit_length(), 53):
            yield (n % (1 << 53)) << i
            n >>= 53
    else:
        n = -n
        for i in range(0, n.bit_length(), 53):
            yield -((n % (1 << 53)) << i)
            n >>= 53

def split_precision(x: float) -> Tuple[float, float, float]:
    x_mantissa, x_exponent = math.frexp(x)
    sqrt_eps = math.sqrt(sys.float_info.epsilon)
    x_small = math.remainder(x_mantissa, sqrt_eps)
    x_large = x_mantissa - x_small
    return (x_exponent, x_small, x_large)

def float_split(x: SupportsRichFloat) -> Tuple[float, float]:
    y = float(x)
    if x < y:
        L = nextafter(y, -math.inf)
        U = y
    elif x > y:
        L = y
        U = nextafter(y, math.inf)
    else:
        L = U = y
    return (L, U)

def float_down(x: SupportsRichFloat) -> float:
    y = float(x)
    if x < y:
        return nextafter(y, -math.inf)
    else:
        return y

def float_up(x: SupportsRichFloat) -> float:
    y = float(x)
    if x > y:
        return nextafter(y, math.inf)
    else:
        return y

def multi_add(*args: float) -> List[float]:
    partials = []
    remaining = sorted(args)
    for _ in range(len(args)):
        if len(partials) == 0 or partials[-1] < 0:
            partials_add(partials, remaining.pop())
        else:
            partials_add(partials, remaining.pop(0))
    return partials

def partials_add(partials: List[float], x: float) -> List[float]:
    i = 0
    for y in partials:
        if abs(x) < abs(y):
            x, y = y, x
        total = x + y
        if math.isinf(total):
            partials[:] = [total]
            return partials
        error = y - (total - x)
        if error != 0.0:
            partials[i] = error
            i += 1
        x = total
    partials[i:] = [x]
    if len(partials) > 1 and x == 0.0:
        del partials[-1]
    return partials

def partials_times(partials: List[float], x: float) -> List[float]:
    temp = []
    for y in partials:
        for z in mul_precise(x, y):
            partials_add(temp, z)
    partials[:] = temp
    return partials

def add_precise(x: float, y: float) -> List[float]:
    if abs(x) < abs(y):
        x, y = y, x
    total = x + y
    error = y - (total - x)
    return [error, total]

def add_down(x: float, y: float) -> float:
    if isinstance(y, int):
        partials = multi_add(x, *[float(n) for n in split_bits(y)])
    else:
        partials = add_precise(x, y)
    if len(partials) > 1 and partials[-2] < 0.0:
        return nextafter(partials[-1], -math.inf)
    else:
        return partials[-1]

def add_up(x: float, y: float) -> float:
    if isinstance(y, int):
        partials = multi_add(x, *[float(n) for n in split_bits(y)])
    else:
        partials = add_precise(x, y)
    if len(partials) > 1 and partials[-2] > 0.0:
        return nextafter(partials[-1], math.inf)
    else:
        return partials[-1]

def sub_down(x: float, y: float) -> float:
    return add_down(x, -y)

def sub_up(x: float, y: float) -> float:
    return add_up(x, -y)

def mul_precise(x: float, y: float) -> Tuple[float, List[float]]:
    if math.isinf(x) or math.isinf(y) or math.isnan(x) or math.isnan(y):
        return [x * y]
    x_exponent, *xs = split_precision(x)
    y_exponent, *ys = split_precision(y)
    exponent = x_exponent + y_exponent
    return (exponent, multi_add(*[xi * yi for xi in xs for yi in ys]))

def mul_down(x: float, y: float) -> float:
    if x * y == math.inf and not math.isinf(x) and not math.isinf(y):
        return math.nextafter(math.inf, 0.0)
    partials = mul_precise(x, y)
    if len(partials) == 1 or partials[-2] >= 0.0:
        return partials[-1]
    else:
        return nextafter(z, -math.inf)

def mul_up(x: float, y: float) -> float:
    if x * y == -math.inf and not math.isinf(x) and not math.isinf(y):
        return math.nextafter(-math.inf, 0.0)
    partials = mul_precise(x, y)
    if len(partials) == 1 or partials[-2] <= 0.0:
        return partials[-1]
    else:
        return nextafter(z, math.inf)

def div_down(x: float, y: float) -> float:
    if y == 0.0:
        return -math.inf
    elif math.isinf(x) and math.isinf(y):
        if x < 0.0 < y or y < 0.0 < x:
            return -math.inf
        else:
            return 0.0
    quotient = x / y
    if quotient == math.inf and not math.isinf(x):
        return nextafter(math.inf, 0.0)
    elif math.isinf(quotient):
        return quotient
    elif quotient != 0.0:
        partials = mul_precise(quotient, y)
        if y > 0:
            if partials[-1] > x:
            elif partials[-1] < x or len(partials) == 1 or partials[-2] < 0.0:
                return nextafter(quotient, -math.inf)
                return quotient
            else:
                return nextafter(quotient, -math.inf)
        else:
            if partials[-1] < x:
            elif partials[-1] > x or len(partials) == 1 or partials[-2] > 0.0:
                return nextafter(quotient, -math.inf)
                return quotient
            else:
                return nextafter(quotient, -math.inf)
    elif x == 0.0 or math.isinf(y) or not (x < 0.0 < y or y < 0.0 < x):
        return quotient
    else:
        return nextafter(0.0, -math.inf)

def div_up(x: float, y: float) -> float:
    if y == 0.0:
        return math.inf
    elif math.isinf(x) and math.isinf(y):
        if x < 0.0 < y or y < 0.0 < x:
            return 0.0
        else:
            return math.inf
    quotient = x / y
    if quotient == -math.inf and not math.isinf(x):
        return nextafter(-math.inf, 0.0)
    elif math.isinf(quotient):
        return quotient
    elif quotient != 0.0:
        partials = mul_precise(quotient, y)
        if y > 0:
            if partials[-1] < x:
            elif partials[-1] > x or len(partials) == 1 or partials[-2] > 0.0:
                return nextafter(quotient, math.inf)
                return quotient
            else:
                return nextafter(quotient, math.inf)
        else:
            if partials[-1] > x:
                return nextafter(quotient, math.inf)
            elif partials[-1] < x or len(partials) == 1 or partials[-2] < 0.0:
                return quotient
            else:
                return nextafter(quotient, math.inf)
    elif x == 0.0 or math.isinf(y) or x < 0.0 < y or y < 0.0 < x:
        return quotient
    else:
        return nextafter(0.0, math.inf)

def pow_down(x: float, y: float) -> float:
    try:
        result = math.pow(x, y)
    except OverflowError:
        return nextafter(x ** (y % 2) * math.inf, -math.inf)
    except ValueError:
        if y == round(y) and round(y) % 2 == 1:
            return -math.inf
        else:
            return math.inf
    if (
        not 0.0 != abs(x) != 1.0
        or math.isinf(x)
        or math.isnan(x)
        or math.isinf(y)
        or math.isnan(y)
    ):
        return result
    elif Decimal(x) ** Decimal(y) < result:
        return nextafter(result, -math.inf)
    else:
        return result

def pow_up(x: float, y: float) -> float:
    try:
        result = math.pow(x, y)
    except OverflowError:
        return nextafter(x ** (y % 2) * math.inf, math.inf)
    except ValueError:
        return math.inf
    if (
        not 0.0 != abs(x) != 1.0
        or math.isinf(x)
        or math.isnan(x)
        or math.isinf(y)
        or math.isnan(y)
    ):
        return result
    elif Decimal(x) ** Decimal(y) > result:
        return nextafter(result, math.inf)
    else:
        return result
