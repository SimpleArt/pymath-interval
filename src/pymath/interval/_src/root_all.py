import math
import sys
from itertools import chain
from math import exp, log, sqrt
from typing import Literal, Optional, TypeVar

if sys.version_info < (3, 9):
    from typing import Callable, Iterator, Tuple
else:
    from builtins import tuple as Tuple
    from collections.abc import Callable, Iterator

T = TypeVar("T")

if sys.version_info < (3, 10):
    if sys.version_info < (3, 9):
        from typing import Iterable
    else:
        from collections.abc import Iterable
    def pairwise(iterable: Iterable[T]) -> Iterator[Tuple[T, T]]:
        """Generates consecutive pairs of elements."""
        iterator = iter(iterable)
        x = next(iterator, None)
        for y in iterator:
            yield (x, y)
            x = y
else:
    from itertools import pairwise

from .fpu_rounding import nextafter
from .interval import Interval

__all__ = ["bisect", "newton"]

def is_between(lo: float, x: float, hi: float, /) -> bool:
    """Checks if `x` is between the `lo` and `hi` arguments."""
    return lo < x < hi or lo > x > hi

def mean(x1: float, x2: float, /) -> float:
    """Returns the arithmetic mean of x1 and x2 without overflowing."""
    return x1 + 0.5 * (x2 - x1) if sign(x1) == sign(x2) else 0.5 * (x1 + x2)

def sign(x: float) -> Literal[-1, 0, 1]:
    """Returns the sign of a real number: -1, 0, or 1."""
    return 1 if x > 0 else -1 if x < 0 else 0

def bisect(
    f: Callable[[Interval], Interval],
    interval: Interval,
    x: Optional[float],
    abs_err: float,
    rel_err: float,
    abs_tol: float,
    rel_tol: float,
    /,
) -> Iterator[Interval]:
    if 0 in f(interval[:nextafter(-math.inf, 0.0)]):
        yield interval[:nextafter(-math.inf, 0.0)]
    leftover = interval[nextafter(math.inf, 0.0):]
    interval = interval[nextafter(-math.inf, 0.0):nextafter(math.inf, 0.0)]
    if x is None:
        intervals = [*interval.sub_intervals]
    else:
        x = float(x)
        intervals = [*interval[:x].sub_intervals, *interval[x:].sub_intervals]
    intervals = [
        interval
        for interval in reversed(intervals)
        if 0 in f(interval)
    ]
    while len(intervals) > 0:
        interval = intervals.pop()
        x1 = interval.minimum
        x2 = interval.maximum
        x = 0.5 * (_utils.sign(x1) + _utils.sign(x2))
        while (
            x1 <= x <= x2
            and not _utils.is_between(x1 / 8, x2, x1 * 8)
            and not (_utils.sign(x1) == _utils.sign(x2) and _utils.is_between(sqrt(abs(x1)), abs(x2), x1 * x1))
            and x2 - x1 > abs_err + rel_err * abs(interval).minimum
        ):
            if abs(x1 - x2) < 16 * (abs_tol + rel_tol * abs(x)):
                x += 0.25 * (abs_err + rel_err * abs(x)) * _utils.sign((x1 - x) + (x2 - x))
            else:
                x += 0.25 * (abs_tol + rel_tol * abs(x)) * _utils.sign((x1 - x) + (x2 - x))
            left = interval[:x]
            right = interval[x:]
            if 0 in f(left):
                if 0 in f(right):
                    intervals.append(right)
                interval = left
            elif 0 in f(right):
                interval = right
            else:
                interval = None
                break
            x1 = interval.minimum
            x2 = interval.maximum
            x = 0.5 * (_utils.sign(x1) + _utils.sign(x2))
            if abs(x) == 0.5:
                x = 0.0
        if interval is None:
            continue
        x_sign = 1 if x > 0 else -1
        x_abs = 1 if abs(x1) > 1 else -1
        while (
            not _utils.is_between(x1 / 8, x2, x1 * 8)
            and not _utils.is_between(sqrt(abs(x1)), abs(x2), x1 * x1)
            and x2 - x1 > abs_err + rel_err * abs(interval).minimum
        ):
            x = x_sign * exp(x_abs * sqrt(log(abs(x1)) * log(abs(x2))))
            if abs(x1 - x2) < 16 * (abs_tol + rel_tol * abs(x)):
                x += 0.25 * (abs_err + rel_err * abs(x)) * _utils.sign((x1 - x) + (x2 - x))
            else:
                x += 0.25 * (abs_tol + rel_tol * abs(x)) * _utils.sign((x1 - x) + (x2 - x))
            left = interval[:x]
            right = interval[x:]
            if 0 in f(left):
                if 0 in f(right):
                    intervals.append(right)
                interval = left
            elif 0 in f(right):
                interval = right
            else:
                interval = None
                break
            x1 = interval.minimum
            x2 = interval.maximum
        if interval is None:
            continue
        while (
            not _utils.is_between(x1 / 8, x2, x1 * 8)
            and x2 - x1 > abs_err + rel_err * abs(interval).minimum
        ):
            x = x_sign * sqrt(abs(x1)) * sqrt(abs(x2))
            if abs(x1 - x2) < 16 * (abs_tol + rel_tol * abs(x)):
                x += 0.25 * (abs_err + rel_err * abs(x)) * _utils.sign((x1 - x) + (x2 - x))
            else:
                x += 0.25 * (abs_tol + rel_tol * abs(x)) * _utils.sign((x1 - x) + (x2 - x))
            left = interval[:x]
            right = interval[x:]
            if 0 in f(left):
                if 0 in f(right):
                    intervals.append(right)
                interval = left
            elif 0 in f(right):
                interval = right
            else:
                interval = None
                break
            x1 = interval.minimum
            x2 = interval.maximum
        if interval is None:
            continue
        while x2 - x1 > abs_err + rel_err * abs(interval).minimum:
            x = x1 + 0.5 * (x2 - x1)
            left = interval[:x]
            right = interval[x:]
            if 0 in f(left):
                if 0 in f(right):
                    intervals.append(right)
                interval = left
            elif 0 in f(right):
                interval = right
            else:
                interval = None
                break
            x1 = interval.minimum
            x2 = interval.maximum
        if interval is not None:
            yield interval
    if 0.0 in f(leftover):
        yield leftover

def newton(
    f: Callable[[Interval], Interval],
    fprime: Callable[[Interval], Interval],
    interval: Interval,
    x: Optional[float],
    abs_err: float,
    rel_err: float,
    abs_tol: float,
    rel_tol: float,
    /,
) -> Iterator[Interval]:
    if 0 in f(interval[:nextafter(-math.inf, 0.0)]):
        yield interval[:nextafter(-math.inf, 0.0)]
    leftover = interval[nextafter(math.inf, 0.0):]
    interval = interval[nextafter(-math.inf, 0.0):nextafter(math.inf, 0.0)]
    if x is None:
        intervals = [*interval.sub_intervals]
    else:
        x = float(x)
        intervals = [*interval[:x].sub_intervals, *interval[x:].sub_intervals]
    intervals = [
        interval
        for interval in reversed(intervals)
        if 0 in f(interval)
    ]
    use_left = True
    while len(intervals) > 0:
        interval = intervals.pop()
        size = interval.size
        if size <= abs_err + rel_err * abs(interval).minimum:
            yield interval
            continue
        if use_left:
            x = interval.minimum
        else:
            x = interval.maximum
        use_left ^= True
        y = f(Interval((x, x)))
        for right, left in pairwise(chain(
            [interval.maximum],
            reversed((interval & (x - y / fprime(interval)))._endpoints),
            [interval.minimum],
        )):
            if left == right:
                continue
            elif left == interval.minimum and right == interval.maximum:
                x = _utils.mean(left, right)
                intervals.extend(
                    sub_interval
                    for sub_interval in [interval[x:], interval[:x]]
                    if 0 in f(sub_interval)
                )
            elif 0 in f(Interval((left, right))):
                interval = Interval((left, right))
                if interval.size <= abs_tol + rel_tol * abs(interval).minimum:
                    abs_current = abs_err
                    rel_current = rel_err
                else:
                    abs_current = abs_tol
                    rel_current = rel_tol
                x = right - 0.25 * (abs_current + rel_current * abs(right))
                x = max(x, _utils.mean(left, right))
                if x < right:
                    if 0 in f(interval[x:]):
                        intervals.append(interval[x:])
                    interval = interval[:x]
                    if interval == interval[()]:
                        continue
                x = left + 0.25 * (abs_current + rel_current * abs(left))
                x = min(x, _utils.mean(left, right))
                iterator = reversed(interval[x:]._endpoints)
                intervals.extend(
                    interval
                    for upper, lower in zip(iterator, iterator)
                    for interval in [Interval((lower, upper))]
                    if 0 in f(interval)
                )
                if x > left:
                    interval = interval[:x]
                    if 0 in f(interval):
                        intervals.append(interval)
    if 0.0 in f(leftover):
        yield leftover
