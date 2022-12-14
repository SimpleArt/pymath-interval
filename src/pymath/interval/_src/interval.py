from __future__ import annotations
import operator
import re
import sys
from decimal import Decimal
from heapq import merge
from math import ceil, floor, inf, isinf, isnan
from typing import Any, Optional, Protocol, SupportsFloat
from typing import SupportsIndex, TypeVar, Union, get_args, overload

if sys.version_info < (3, 9):
    from typing import Iterator, Pattern, Tuple, Type
else:
    from builtins import tuple as Tuple, type as Type
    from collections.abc import Iterator
    from re import Pattern

from . import fpu_rounding as fpur
from .typing import RealLike, SupportsRichFloat

__all__ = ["Interval", "interval"]

Self = TypeVar("Self", bound="Interval")
SupportsSelf = TypeVar("SupportsSelf", bound="SupportsInterval")

NOT_REAL = "could not interpret {} as a real value"

FSTRING_FORMATTER: Pattern = re.compile(
    "(?P<fill>.*?)"
    "(?P<align>[<>=^]?)"
    "(?P<sign>[+ -]?)"
    "(?P<alternate>[#]?)"
    "(?P<width>[0-9]*)"
    "(?P<group>[_,]?)"
    "(?P<precision>(?:[.][0-9]+)?)"
    "(?P<dtype>[bcdeEfFgGnosxX%]?)"
)


class SupportsInterval(Protocol):

    def __or__(self: SupportsSelf, other: Interval) -> SupportsSelf: ...

    def __ror__(
        self: SupportsSelf,
        other: Union[Interval, RealLike, SupportsSelf],
    ) -> SupportsSelf: ...


class Interval:
    _endpoints: tuple[float, ...]

    __slots__ = ("_endpoints",)

    def __init__(self: Self, /, *args: Tuple[RealLike, RealLike]) -> None:
        for arg in args:
            if not isinstance(arg, tuple):
                raise TypeError(f"interval(...) expects tuples for arguments, got {arg!r}")
            elif len(arg) != 2:
                raise ValueError(f"interval(...) expects (lower, upper) for arguments, got {len(arg)!r} arguments")
            for x in arg:
                if not isinstance(x, get_args(RealLike)):
                    raise TypeError(NOT_REAL.format(repr(x)))
                elif isinstance(x, Decimal) and x.is_nan():
                    raise TypeError(NOT_REAL.format(repr(x)))
                elif isinstance(x, SupportsFloat) and isnan(float(x)):
                    raise ValueError(NOT_REAL.format(repr(float(x))))
        intervals = []
        for lower, upper in args:
            if isinstance(lower, SupportsIndex):
                lower = operator.index(lower)
            lower = fpur.float_down(lower)
            if isinstance(upper, SupportsIndex):
                upper = operator.index(upper)
            upper = fpur.float_up(upper)
            if lower <= upper:
                intervals.append((lower, upper))
        intervals.sort()
        if len(intervals) == 0:
            self._endpoints = ()
            return
        endpoints = [intervals[0][0]]
        upper = intervals[0][1]
        for L, U in intervals:
            if L > upper:
                endpoints.append(upper)
                endpoints.append(L)
                upper = U
            elif upper < U:
                upper = U
        endpoints.append(upper)
        self._endpoints = (*endpoints,)

    def __abs__(self: Self, /) -> Self:
        return -self[:0] | self[0:]

    def __add__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        iterator = iter(self._endpoints)
        if isinstance(other, Interval) and type(self).__add__ is type(other).__add__:
            return type(self)(*[
                (fpur.add_down(x_lower, y_lower), fpur.add_up(x_upper, y_upper))
                for x_lower, x_upper in zip(iterator, iterator)
                for y_lower, y_upper in zip(*[iter(other._endpoints)] * 2)
            ])
        elif isinstance(other, Decimal):
            return type(self)(*[
                (fpur.float_down(Decimal(lower) + other), fpur.float_up(Decimal(upper) + other))
                for lower, upper in zip(iterator, iterator)
            ])
        elif isinstance(other, SupportsIndex):
            other = operator.index(other)
            return type(self)(*[
                (fpur.add_down(lower, other), fpur.add_up(upper, other))
                for lower, upper in zip(iterator, iterator)
            ])
        elif isinstance(other, SupportsFloat):
            return self + Interval(fpur.float_split(other))
        else:
            return NotImplemented

    def __and__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        iterator = iter(self._endpoints)
        if isinstance(other, Interval) and type(self).__and__ is type(other).__and__:
            return type(self)(*[
                (max(x_lower, y_lower), min(x_upper, y_upper))
                for x_lower, x_upper in zip(iterator, iterator)
                for y_lower, y_upper in zip(*[iter(other._endpoints)] * 2)
            ])
        elif isinstance(other, get_args(RealLike)):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return self & Interval(fpur.float_split(other))
        else:
            return NotImplemented

    def __as_interval__(self: Self) -> Interval:
        return self

    @overload
    def __call__(
        self: Self,
        /,
        *args: Union[Interval, RealLike],
    ) -> Self: ...

    @overload
    def __call__(
        self: Self,
        arg: SupportsSelf,
        /,
        *args: Union[Interval, RealLike, SupportsSelf],
    ) -> SupportsSelf: ...

    def __call__(self, /, *args):
        result = self[()]
        for arg in args:
            result = arg | result
        return result

    def __contains__(self: Self, other: Any, /) -> bool:
        if isinstance(other, SupportsIndex):
            other = operator.index(other)
        elif isinstance(other, SupportsFloat):
            return all(
                any(
                    x.minimum <= o <= x.maximum
                    for x in self.sub_intervals
                )
                for o in {*fpur.float_split(other)}
            )
        elif not isinstance(other, Decimal):
            return False
        return any(
            x.minimum <= other <= x.maximum
            for x in self.sub_intervals
        )

    @classmethod
    def __dist__(cls: Type[Self], p: List[Interval], q: List[Interval]) -> Self:
        return NotImplemented

    def __eq__(self: Self, other: Any, /) -> bool:
        if isinstance(other, Interval) and type(self).__eq__ is type(other).__eq__:
            return self._endpoints == other._endpoints
        elif isinstance(other, (Decimal, SupportsIndex)):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return self._endpoints == (other, other)
        elif isinstance(other, SupportsFloat):
            return self._endpoints == fpur.float_split(other)
        else:
            return NotImplemented

    def __format__(self: Self, specifier: str, /) -> str:
        match = FSTRING_FORMATTER.fullmatch(specifier)
        if match is None:
            # Invalid specifier, raise an error.
            f"{0.0:{specifier}}"
            raise TypeError(f"invalid specifier, got {specifier!r}")
        fill, align, sign, alternate, width, group, precision, dtype = match.groups()
        if len(self._endpoints) == 0:
            return "interval()"
        if dtype == "":
            dtype = "g"
            if precision == "":
                precision = ".17"
        specifier = f"{fill}{align}{sign}{alternate}{width}{group}{precision}{dtype}"
        iterator = iter(self._endpoints)
        if self.size == 0:
            points = ", ".join([f"{U:{specifier}}" for _, U in zip(iterator, iterator)])
            return f"interval({points})"
        bounds = ", ".join([
            ":".join([
                f"{-0.0:{specifier}}" if lower == 0.0 else "" if isinf(lower) and lower < 0.0 else f"{lower:{specifier}}",
                f"{0.0:{specifier}}" if upper == 0.0 else "" if isinf(upper) and upper > 0.0 else f"{upper:{specifier}}",
            ])
            for lower, upper in zip(iterator, iterator)
        ])
        return f"interval[{bounds}]"

    @classmethod
    def __fsum__(cls: Type[Self], intervals: List[Interval]) -> Self:
        return NotImplemented

    def __getitem__(self: Self, args: Union[slice, Tuple[slice, ...]], /) -> Interval:
        if isinstance(args, slice):
            args = (args,)
        elif not isinstance(args, tuple):
            raise TypeError(f"interval[...] expects 0 or more slices, got {args!r}")
        for arg in args:
            if not isinstance(arg, slice):
                raise TypeError(f"interval[...] expects slices, got {arg!r}")
            elif arg.step is not None:
                raise TypeError(f"interval[...] expects [lower:upper] for arguments, got a step argument")
            elif arg.start is not None and not isinstance(arg.start, (Decimal, SupportsIndex, SupportsFloat)):
                raise TypeError(NOT_REAL.format(repr(arg.start)))
            elif arg.start is not None and isinstance(arg.start, Decimal) and arg.start.is_nan():
                raise TypeError(NOT_REAL.format(repr(arg.start)))
            elif arg.start is not None and isinstance(arg.start, SupportsFloat) and isnan(float(arg.start)):
                raise TypeError(NOT_REAL.format(repr(float(arg.start))))
            elif arg.stop is not None and not isinstance(arg.stop, (Decimal, SupportsIndex, SupportsFloat)):
                raise TypeError(NOT_REAL.format(repr(arg.stop)))
            elif arg.stop is not None and isinstance(arg.stop, Decimal) and arg.stop.is_nan():
                raise TypeError(NOT_REAL.format(repr(arg.stop)))
            elif arg.stop is not None and isinstance(arg.stop, SupportsFloat) and isnan(float(arg.stop)):
                raise TypeError(NOT_REAL.format(repr(float(arg.stop))))
        intervals = []
        for arg in args:
            if arg.start is None:
                L = -inf
            elif isinstance(arg.start, SupportsIndex):
                L = fpur.float_down(operator.index(arg.start))
            else:
                L = fpur.float_down(arg.start)
            if arg.stop is None:
                U = inf
            elif isinstance(arg.stop, SupportsIndex):
                U = fpur.float_up(operator.index(arg.stop))
            else:
                U = fpur.float_up(arg.stop)
            intervals.append((L, U))
        return self & type(self)(*intervals)

    def __hash__(self: Self, /) -> int:
        return hash(self._endpoints)

    def __invert__(self: Self, /) -> Self:
        iterator = iter([-inf, *self._endpoints, inf])
        return type(self)(*zip(iterator, iterator))

    def __mul__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval) and type(self).__mul__ is type(other).__mul__:
            intervals = []
            if 0 in self and len(other._endpoints) > 0 or 0 in other and len(self._endpoints) > 0:
                if any(isinf(x) for x in (self.minimum, self.maximum, other.minimum, other.maximum)):
                    return interval
                intervals.append((0, 0))
            for x in self[0:].sub_intervals:
                if x.maximum == 0:
                    continue
                for y in other[0:].sub_intervals:
                    if y.maximum == 0:
                        continue
                    try:
                        start = fpur.mul_down(x.minimum, y.minimum)
                    except OverflowError:
                        start = inf
                    try:
                        stop = fpur.mul_up(x.maximum, y.maximum)
                    except OverflowError:
                        intervals.append((start, inf))
                    else:
                        intervals.append((start, stop))
                for y in other[:0].sub_intervals:
                    if y.minimum == 0:
                        continue
                    try:
                        stop = fpur.mul_up(x.minimum, y.maximum)
                    except OverflowError:
                        stop = -inf
                    try:
                        start = fpur.mul_down(x.maximum, y.minimum)
                    except OverflowError:
                        intervals.append((-inf, stop))
                    else:
                        intervals.append((start, stop))
            for x in self[:0].sub_intervals:
                if x.minimum == 0:
                    continue
                for y in other[:0].sub_intervals:
                    if y.minimum == 0:
                        continue
                    try:
                        start = fpur.mul_down(x.maximum, y.maximum)
                    except OverflowError:
                        start = inf
                    try:
                        stop = fpur.mul_up(x.minimum, y.minimum)
                    except OverflowError:
                        intervals.append((start, inf))
                    else:
                        intervals.append((start, stop))
                for y in other[0:].sub_intervals:
                    if y.maximum == 0:
                        continue
                    try:
                        stop = fpur.mul_up(x.maximum, y.minimum)
                    except OverflowError:
                        stop = -inf
                    try:
                        start = fpur.mul_down(x.minimum, y.maximum)
                    except OverflowError:
                        intervals.append((-inf, stop))
                    else:
                        intervals.append((start, stop))
            return type(self)(*intervals)
        elif isinstance(other, get_args(RealLike)):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return self * Interval(fpur.float_split(other))
        else:
            return NotImplemented

    def __neg__(self: Self, /) -> Self:
        iterator = reversed(self._endpoints)
        return type(self)(*[(-upper, -lower) for upper, lower in zip(iterator, iterator)])

    def __or__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval) and type(self).__or__ is type(other).__or__:
            return type(self)(
                *[(x.minimum, x.maximum) for x in self.sub_intervals],
                *[(x.minimum, x.maximum) for x in other.sub_intervals],
            )
        elif isinstance(other, get_args(RealLike)):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return type(self)(*[(x.minimum, x.maximum) for x in self.sub_intervals], fpur.float_split(other))
        else:
            return NotImplemented

    def __pos__(self: Self, /) -> Self:
        return self

    def __pow__(self: Self, other: Union[Interval, RealLike], modulo: None = None, /) -> Interval:
        if modulo is not None:
            return NotImplemented
        elif isinstance(other, Interval) and type(self).__pow__ is type(other).__pow__:
            intervals = []
            for x in self[0:].sub_intervals:
                for y in other[0:].sub_intervals:
                    start = fpur.pow_down(x.minimum, y.minimum)
                    stop = fpur.pow_up(x.maximum, y.maximum)
                    intervals.append((start, stop))
                for y in other[:0].sub_intervals:
                    if y.minimum < 0.0 == x.minimum:
                        start = stop = inf
                    else:
                        start = fpur.pow_down(x.maximum, y.maximum)
                        stop = fpur.pow_up(x.minimum, y.minimum)
                    intervals.append((start, stop))
            return type(self)(*intervals)
        elif isinstance(other, SupportsIndex):
            other = operator.index(other)
            intervals = []
            iterator = iter(self._endpoints)
            if other == 0:
                for _ in iterator:
                    intervals.append((1.0, 1.0))
                    break
            elif other > 0 and other % 2 == 0:
                for lower, upper in zip(iterator, iterator):
                    if upper < 0:
                        intervals.append((fpur.pow_down(upper, other), fpur.pow_up(lower, other)))
                    elif lower > 0:
                        intervals.append((fpur.pow_down(lower, other), fpur.pow_up(upper, other)))
                    else:
                        intervals.append((0.0, fpur.pow_up(max(lower, upper, key=abs), other)))
            elif other > 0:
                intervals = [
                    (fpur.pow_down(lower, other), fpur.pow_up(upper, other))
                    for lower, upper in zip(iterator, iterator)
                ]
            elif other % 2 == 0:
                for lower, upper in zip(iterator, iterator):
                    if upper < 0:
                        intervals.append((fpur.pow_down(lower, other), fpur.pow_up(upper, other)))
                    elif lower > 0:
                        intervals.append((fpur.pow_down(upper, other), fpur.pow_up(lower, other)))
                    elif lower == upper:
                        intervals.append((inf, inf))
                    else:
                        intervals.append((fpur.pow_down(max(lower, upper, key=abs), other), inf))
            else:
                for lower, upper in zip(iterator, iterator):
                    if not lower <= 0 <= upper:
                        intervals.append((fpur.pow_down(upper, other), fpur.pow_up(lower, other)))
                        continue
                    elif lower == 0 == upper:
                        intervals.append((-inf, -inf))
                        intervals.append((inf, inf))
                        continue
                    if lower < 0:
                        intervals.append((-inf, fpur.pow_up(lower, other)))
                    if upper > 0:
                        intervals.append((fpur.pow_down(upper, other), inf))
            return type(self)(*intervals)
        elif isinstance(other, SupportsFloat):
            other = Interval(fpur.float_split(other))
            if other.minimum > 0:
                iterator = iter(self[0:]._endpoints)
                intervals = [
                    (fpur.pow_down(lower, other.minimum), fpur.pow_up(upper, other.maximum))
                    for lower, upper in zip(iterator, iterator)
                ]
                return type(self)(*intervals)
            else:
                iterator = iter(self[0:]._endpoints)
                intervals = [
                    (fpur.pow_down(upper, other.maximum), fpur.pow_up(lower, other.minimum))
                    for lower, upper in zip(iterator, iterator)
                ]
                return type(self)(*intervals)
        else:
            return NotImplemented

    def __radd__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval):
            return other.__as_interval__() + self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            return self + other
        else:
            return NotImplemented

    def __rand__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval):
            return other.__as_interval__() & self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            return self & other
        else:
            return NotImplemented

    def __repr__(self: Self, /) -> str:
        return f"{self}"

    def __rmul__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval):
            return other.__as_interval__() * self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            return self * other
        else:
            return NotImplemented

    def __ror__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval):
            return other.__as_interval__() | self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            return self | other
        else:
            return NotImplemented

    def __rpow__(self: Self, other: Union[Interval, RealLike], modulo: None = None, /) -> Interval:
        if modulo is not None:
            return NotImplemented
        elif isinstance(other, Interval):
            return other.__as_interval__() ** self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return Interval(fpur.float_split(other)) ** self
        else:
            return NotImplemented

    def __rsub__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, (Interval, *get_args(RealLike))):
            return -self + other
        else:
            return NotImplemented

    def __rtruediv__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval):
            return other.__as_interval__() / self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return Interval(fpur.float_split(other)) / self
        else:
            return NotImplemented

    def __rxor__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval):
            return other.__as_interval__() ^ self.__as_interval__()
        elif isinstance(other, get_args(RealLike)):
            return self ^ other
        else:
            return NotImplemented

    def __sub__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        iterator = iter(self._endpoints)
        if isinstance(other, Interval) and type(self).__sub__ is type(other).__sub__:
            return type(self)(*[
                (fpur.sub_down(x_lower, y_upper), fpur.sub_up(x_upper, y_lower))
                for x_lower, x_upper in zip(iterator, iterator)
                for y_lower, y_upper in zip(*[iter(other._endpoints)] * 2)
            ])
        elif isinstance(other, Interval):
            return -other + self
        elif isinstance(other, Decimal):
            return type(self)(*[
                (fpur.float_down(Decimal(lower) - other), fpur.float_up(Decimal(upper) - other))
                for lower, upper in zip(iterator, iterator)
            ])
        elif isinstance(other, SupportsIndex):
            other = operator.index(other)
            return type(self)(*[
                (fpur.sub_down(lower, other), fpur.sub_up(upper, other))
                for lower, upper in zip(iterator, iterator)
            ])
        elif isinstance(other, SupportsFloat):
            return self - Interval(fpur.float_split(other))
        else:
            return NotImplemented

    def __truediv__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval) and type(self).__truediv__ is type(other).__truediv__:
            intervals = []
            if not len(self._endpoints) != 0 != len(other._endpoints):
                return type(self)()
            if 0 in self:
                if 0 in other:
                    return type(self)((-inf, inf))
                intervals.append((0, 0))
            elif 0 in other:
                intervals.append((-inf, -inf))
                intervals.append((inf, inf))
            for x in self[0:].sub_intervals:
                if x.maximum == 0:
                    continue
                for y in other[0:].sub_intervals:
                    if y.maximum == 0:
                        continue
                    try:
                        start = fpur.div_down(x.minimum, y.maximum)
                    except (OverflowError, ZeroDivisionError):
                        start = inf
                    try:
                        stop = fpur.div_up(x.maximum, y.minimum)
                    except (OverflowError, ZeroDivisionError):
                        intervals.append((start, inf))
                    else:
                        intervals.append((start, stop))
                for y in other[:0].sub_intervals:
                    if y.minimum == 0:
                        continue
                    try:
                        stop = fpur.div_up(x.minimum, y.minimum)
                    except (OverflowError, ZeroDivisionError):
                        stop = -inf
                    try:
                        start = fpur.div_down(x.maximum, y.maximum)
                    except (OverflowError, ZeroDivisionError):
                        intervals.append((-inf, stop))
                    else:
                        intervals.append((start, stop))
            for x in self[:0].sub_intervals:
                if x.minimum == 0:
                    continue
                for y in other[:0].sub_intervals:
                    if y.minimum == 0:
                        continue
                    try:
                        start = fpur.div_down(x.maximum, y.minimum)
                    except (OverflowError, ZeroDivisionError):
                        start = inf
                    try:
                        stop = fpur.div_up(x.minimum, y.maximum)
                    except (OverflowError, ZeroDivisionError):
                        intervals.append((start, inf))
                    else:
                        intervals.append((start, stop))
                for y in other[0:].sub_intervals:
                    if y.maximum == 0:
                        continue
                    try:
                        stop = fpur.div_up(x.maximum, y.maximum)
                    except (OverflowError, ZeroDivisionError):
                        stop = -inf
                    try:
                        start = fpur.div_down(x.minimum, y.minimum)
                    except (OverflowError, ZeroDivisionError):
                        intervals.append((-inf, stop))
                    else:
                        intervals.append((start, stop))
            return type(self)(*intervals)
        elif isinstance(other, get_args(RealLike)):
            if isinstance(other, SupportsFloat):
                other = operator.index(other)
            other = Interval(fpur.float_split(other))
            if other == 0:
                if 0 in self:
                    return type(self)((-inf, inf))
                elif self == type(self)():
                    return self
                else:
                    return type(self)((-inf, -inf), (inf, inf))
            elif other.minimum >= 0:
                iterator = iter(self._endpoints)
                return type(self)(*[
                    (fpur.div_down(lower, L), fpur.div_up(upper, U))
                    for lower, upper in zip(iterator, iterator)
                    for L in [other.minimum if lower < 0.0 else other.maximum]
                    for U in [other.maximum if upper < 0.0 else other.minimum]
                ])
            elif other.maximum <= 0:
                iterator = reversed(self._endpoints)
                return type(self)(*[
                    (fpur.div_down(upper, L), fpur.div_up(lower, U))
                    for upper, lower in zip(iterator, iterator)
                    for L in [other.minimum if upper < 0.0 else other.maximum]
                    for U in [other.maximum if lower < 0.0 else other.minimum]
                ])
            else:
                assert False, "float_split should never produce values separated by 0"
        else:
            return NotImplemented

    def __xor__(self: Self, other: Union[Interval, RealLike], /) -> Interval:
        if isinstance(other, Interval) and type(self).__xor__ is type(other).__xor__:
            iterator = merge(self._endpoints, other._endpoints)
            return type(self)(*zip(iterator, iterator))
        elif isinstance(other, RealLike):
            if isinstance(other, SupportsIndex):
                other = operator.index(other)
            return self ^ type(self)(fpur.float_split(other))
        else:
            return NotImplemented

    def arccos(self: Self) -> Interval:
        from .imath import acos
        return acos(self)

    def arccosh(self: Self) -> Interval:
        from .imath import acosh
        return acosh(self)

    def arcsin(self: Self) -> Interval:
        from .imath import asin
        return asin(self)

    def arcsinh(self: Self) -> Interval:
        from .imath import asinh
        return asinh(self)

    def arctan(self: Self) -> Interval:
        from .imath import atan
        return atan(self)

    def arctan2(self: Self, other: Union[Interval, RealLike]) -> Interval:
        from .imath import atan2
        return atan2(self, other)

    def arctanh(self: Self) -> Interval:
        from .imath import atanh
        return atanh(self)

    def cos(self: Self) -> Interval:
        from .imath import cos
        return cos(self)

    def cosh(self: Self) -> Interval:
        from .imath import cosh
        return cosh(self)

    def degrees(self: Self, /) -> Interval:
        from .imath import degrees
        return degrees(self)

    rad2deg = degrees

    def exp(self: Self, /) -> Interval:
        from .imath import exp
        return exp(self)

    def exp2(self: Self, /) -> Interval:
        return 2 ** self

    def expm1(self: Self, /) -> Interval:
        from .imath import expm1
        return expm1(self)

    def hypot(self: Self, other: Union[Interval, RealLike]) -> Interval:
        from .imath import hypot
        return hypot(self, other)

    def log(self: Self, base: Union[Interval, RealLike], /) -> Interval:
        from .imath import log
        return log(self, base)

    def log10(self: Self, /) -> Interval:
        from .imath import log10
        return log10(self)

    def log1p(self: Self, /) -> Interval:
        from .imath import log1p
        return log1p(self)

    def log2(self: Self, /) -> Interval:
        from .imath import log2
        return log2(self)

    def radians(self: Self, /) -> Interval:
        from .imath import radians
        return radians(self)

    deg2rad = radians

    def reciprocal(self: Self, /) -> Interval:
        return 1 / self

    def sin(self: Self, /) -> Interval:
        from .imath import sin
        return sin(self)

    def sinh(self: Self, /) -> Interval:
        from .imath import sinh
        return sinh(self)

    def sqrt(self: Self, /) -> Interval:
        from .imath import sqrt
        return sqrt(self)

    def square(self: Self, /) -> Interval:
        return self ** 2

    def tan(self: Self, /) -> Interval:
        from .imath import tan
        return tan(self)

    def tanh(self: Self, /) -> Interval:
        from .imath import tanh
        return tanh(self)

    @property
    def maximum(self: Self, /) -> float:
        if len(self._endpoints) == 0:
            raise ValueError(f"an empty interval has no maximum")
        else:
            return self._endpoints[-1]

    @property
    def minimum(self: Self, /) -> float:
        if len(self._endpoints) == 0:
            raise ValueError(f"an empty interval has no minimum")
        else:
            return self._endpoints[0]

    @property
    def size(self: Self, /) -> float:
        return sum(
            interval.maximum - interval.minimum
            for interval in self.sub_intervals
            if interval.minimum != interval.maximum
        )

    @property
    def sub_intervals(self: Self, /) -> Iterator[Self]:
        iterator = iter(self._endpoints)
        return map(type(self), zip(iterator, iterator))


interval = Interval((-inf, inf))
