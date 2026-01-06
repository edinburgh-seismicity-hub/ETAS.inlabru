# Generate a set of time bins for a specific event.

Generate a set of time bins for a specific event.

## Usage

``` r
time_grid(data.point, coef.t, delta.t, N.exp., T1., T2.)
```

## Arguments

- data.point:

  Point for which the binning is calculated, `list` with elements time
  (`ts, scalar`), event index (`idx.p, scalar`). Names are mandatory and
  should not be changed

- coef.t:

  TimeBinning parameter: look
  [`breaks_exp()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/breaks_exp.md)

- delta.t:

  TimeBinning parameter: look
  [`breaks_exp()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/breaks_exp.md)

- N.exp.:

  TimeBinning parameter: look
  [`breaks_exp()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/breaks_exp.md)

- T1.:

  Start of the temporal domain, `scalar`

- T2.:

  End of the temporal domain `scalar`.

## Value

A `data.frame` with as many rows as the number of bins and fixed number
of columns. The columns are

- `t.start` : starting point of the bin (minimum = `T1.`).

- `t.end` : end point of the bin. (maximum = `T2.`).

- `t.bin.name` : unique bin identifier.

- `t.ref_layer` : bin identifier for calculations

- `ts` : time provided in `data.point`

- `idx.p` : identifier provided in `data.point`

The bins are only between `T1.` and `T2.` or containing `T1.`

## Examples

``` r
## EXAMPLE 1
event <- list(ts = 0, idx.p = 1)
time_grid(data.point = event, coef.t = 1, delta.t = 0.1, N.exp. = 8, T1. = 1, T2. = 20)
#>   t.start t.end t.bin.name t.ref_layer ts idx.p
#> 1     1.0   1.6    0.8-1.6   between-1  0     1
#> 2     1.6   3.2    1.6-3.2           6  0     1
#> 3     3.2   6.4    3.2-6.4           7  0     1
#> 4     6.4  12.8   6.4-12.8           8  0     1
#> 5    12.8  20.0    12.8-20      last-1  0     1
```
