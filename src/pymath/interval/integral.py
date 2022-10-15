import sys
from functools import lru_cache
from math import fsum, inf, isinf, isnan
from typing import Optional, SupportsFloat

if sys.version_info < (3, 9):
    from typing import Callable, Iterator, List, Tuple
else:
    from builtins import list as List, tuple as Tuple
    from collections.abc import Callable, Iterator

from . import imath
from ._src.fpu_rounding import nextafter
from ._src.interval import Interval

def _fmean(x: float, y: float) -> float:
    assert x <= y
    if x <= 0 <= y:
        return 0.5 * (x + y)
    else:
        return x + 0.5 * (y - x)

def _imean(interval: Interval) -> float:
    return _fmean(*interval._endpoints)

def _split(interval: Interval) -> Tuple[Interval, Interval]:
    assert len(interval._endpoints) == 2
    assert interval.minimum < interval.maximum
    if isinf(interval.minimum):
        if interval.maximum > 0:
            return interval[:0], interval[0:]
        else:
            return interval[:interval.maximum * 2 - 1], interval[interval.maximum * 2 - 1:]
    elif isinf(interval.maximum):
        if interval.minimum < 0:
            return interval[:0], interval[0:]
        else:
            return interval[:interval.minimum * 2 + 1], interval[interval.minimum * 2 + 1:]
    else:
        midpoint = _imean(interval)
        return interval[:midpoint], interval[midpoint:]

def partition(
    integral: Callable[[Interval], Interval],
    bounds: Interval,
    error: float = 0.1,
) -> List[Interval]:
    if not callable(integral):
        raise TypeError(f"expected callable integral, got {integral!r}")
    elif not isinstance(bounds, Interval):
        raise TypeError(f"expected interval bounds, got {bounds!r}")
    elif not isinstance(error, SupportsFloat):
        raise TypeError(f"could not interpret error as a real number, got {error!r}")
    elif type(bounds) is not Interval:
        bounds = type(bounds).__as_interval__(bounds)
    error = float(error)
    if isnan(error):
        raise ValueError(f"error cannot be nan")
    partitions = {
        interval: integral(interval)
        for interval in bounds.sub_intervals
        if len(interval._endpoints) > 0
    }
    extremes = {}
    while True:
        if fsum(value.size for value in partitions.values()) < error:
            result = [*partitions, *extremes]
            result.sort(key=lambda p: p.minimum)
            return result
        split_size = max(p.size for p in partitions.values())
        split_size = 0.125 * min(split_size, error)
        split_partitions = {
            k: v
            for k, v in partitions.items()
            if v.size > split_size
        }
        while split_partitions:
            k, v = split_partitions.popitem()
            p = [*k._endpoints]
            if nextafter(*p) == p[1]:
                extremes[k] = partitions.pop(k)
                continue
            elif isinf(p[0]) or isinf(p[1]):
                left, right = _split(k)
                p = [*left._endpoints, right.maximum]
            else:
                if isinf(v.size):
                    n = round((v / 2).size / (split_size / 2)).bit_length()
                else:
                    n = round(v.size / split_size).bit_length()
                for _ in range(n // 3 + 1):
                    p.extend(_fmean(p[i - 1], p[i]) for i in range(1, len(p)))
                    p.sort()
            del partitions[k]
            for i in range(1, len(p)):
                k = Interval((p[i - 1], p[i]))
                v = partitions[k] = integral(k)
                if v.size > split_size:
                    split_partitions[k] = v

def integrate(
    integral: Callable[[Interval], Interval],
    bounds: Interval,
    f: Optional[Callable[[Interval], Interval]] = None,
    error: float = 0.1,
) -> Tuple[Interval, float]:
    if not callable(integral):
        raise TypeError(f"expected callable integral, got {integral!r}")
    elif not isinstance(bounds, Interval):
        raise TypeError(f"expected interval bounds, got {bounds!r}")
    elif f is not None and not callable(f):
        raise TypeError(f"expected None or callable f, got {f!r}")
    elif not isinstance(error, SupportsFloat):
        raise TypeError(f"could not interpret error as a real number, got {error!r}")
    elif type(bounds) is not Interval:
        bounds = type(bounds).__as_interval__(bounds)
    error = float(error)
    if isnan(error):
        raise ValueError(f"error cannot be nan")
    partitions = {
        p: integral(p)
        for p in partition(integral, bounds, error)
    }
    if f is None:
        return imath.fsum(partitions.values()), fsum(map(_imean, partitions.values()))
    else:
        return imath.fsum(partitions.values()), fsum(
            _imean(v)
                if
            isinf(k.size)
            or sum(1 for _ in v.sub_intervals) > 2
            or isinf(v.minimum)
            or isinf(v.maximum)
                else
            _imean(f(k & k.minimum) + 4 * f(k & _imean(k)) + f(k & k.maximum)) / 6
            for k, v in partitions.items()
        )
