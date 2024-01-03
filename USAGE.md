# Usage

## For basic interval arithmetic

```
>>> from pymath.interval import interval

# To intialize an interval
>>> interval(1,5)
interval(1,5)

# To add two intervals
>>> interval(1,5) + interval(7,10)
interval(8, 11, 12, 15)
```
