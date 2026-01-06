# Integral of Omori's law

Function to compute the integral of Omori's law efficiently

## Usage

``` r
compute_grid(param., list.input_)
```

## Arguments

- param.:

  ETAS parameters vector (\\\mu, K, \alpha, c, p\\), only \\c, p\\ are
  used.

- list.input\_:

  `list` containing information to calculate the integrals efficiently.
  The list is created inside the `Temporal.ETAS` function Two elements
  are used

  - `time.sel` selection of rows of the output of `time_grid` with
    unique `t.ref_layer` value, `data.frame`.

  - `Imapping` mapper between the unique names provided in `time.sel`
    and original rows of the output of
    [`time_grid()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/time_grid.md),
    `vector`.

## Value

`vector` with same length as `list.input_$Imapping` with the integral of
Omori's law for each bin.
