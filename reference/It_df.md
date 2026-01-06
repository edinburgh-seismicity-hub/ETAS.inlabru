# Function to calculate the integral of Omori's law

Function to calculate the integral of Omori's law

## Usage

``` r
It_df(param_, time.df)
```

## Arguments

- param\_:

  ETAS parameters vector (\\\mu, K, \alpha, c, p\\), only \\c, p\\ are
  used.

- time.df:

  output of the function
  [`time_grid()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/time_grid.md)

## Value

`vector` of integral values between each bin provided in `time.df`
