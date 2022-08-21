from functools import lru_cache
from math import fsum, inf, isinf, isnan
from typing import Callable, Iterator, Optional, SupportsFloat, Tuple

from . import imath
from ._src.interval import Interval

def _mean(interval: Interval) -> float:
    assert len(interval._endpoints) == 2
    if 0 in interval:
        return 0.5 * (interval.maximum + interval.minimum)
    else:
        return interval.minimum + 0.5 * (interval.maximum - interval.minimum)

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
        midpoint = _mean(interval)
        return interval[:midpoint], interval[midpoint:]

def partition(
    integral: Callable[[Interval], Interval],
    bounds: Interval,
    error: float = 0.1,
) -> Iterator[Interval]:
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
        split_partitions = [
            k
            for k, v in partitions.items()
            if v.size > split_size
        ]
        while split_partitions:
            p = split_partitions.pop()
            left, right = _split(p)
            if left != p != right:
                del partitions[p]
                for p in (left, right):
                    partitions[p] = integral(p)
                    if 8 * partitions[p].size > split_size:
                        split_partitions.append(p)
            else:
                extremes[p] = partitions.pop(p)

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
        return imath.fsum(partitions.values()), fsum(map(_mean, partitions.values()))
    else:
        return imath.fsum(partitions.values()), fsum(
            _mean(v)
                if
            isinf(k.size)
            or sum(1 for _ in v.sub_intervals) > 2
            or isinf(v.minimum)
            or isinf(v.maximum)
                else
            (f(v & k.minimum) + 4 * f(v & _mean(k)) + f(v & k.maximum)) / 6
            for k, v in partitions.items()
        )
