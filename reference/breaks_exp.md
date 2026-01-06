# Find breaks point for 1D grid

`breaks_exp` return the breaks points of a one dimensional grid
depending on three parameters, see details

## Usage

``` r
breaks_exp(start.grid, end.grid, coef.t = 2, delta.t, N.exp. = 10)
```

## Arguments

- start.grid:

  Starting point of the grid, `scalar`.

- end.grid:

  End point of the grid, `scalar`.

- coef.t:

  TimeBinning parameter: \\\delta \> 0\\ determines the relative length
  of subsequent intervals, `scalar`.

- delta.t:

  TimeBinning parameter: \\\Delta \> 0\\ determines the length of the
  intervals, `scalar`.

- N.exp.:

  TimeBinning parameter: \\n\_{max} \>0\\ determines the maximum number
  of intervals, `scalar`

## Value

`vector` containing the grid points

## Details

The grid is calculated as follows \$\$t, t + \Delta, t + \Delta(1 +
\delta), t + \Delta(1 + \delta)^2,...., T\$\$ where \\t\\ is the
`start.grid` argument, \\T\\ is the `end.grid` argument, and
\\n\_{max}\\ is the maximum value of the exponent

## Examples

``` r
breaks_exp(
  start.grid = 1, end.grid = 100, coef.t = 1, delta.t = 1, N.exp. = 3
)
#> [1]   1   2   3   5   9 100
breaks_exp(
  start.grid = 1, end.grid = 100, coef.t = 1, delta.t = 1, N.exp. = 10
)
#> [1]   1   2   3   5   9  17  33  65 100
```
