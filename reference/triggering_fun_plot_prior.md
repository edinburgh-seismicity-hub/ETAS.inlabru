# Triggering function plot from prior samples

Function to plot the ETAS triggering function corresponding to different
prior samples

## Usage

``` r
triggering_fun_plot_prior(
  input.list,
  magnitude = 4,
  n.samp = 10,
  t.end = 1,
  n.breaks = 100
)
```

## Arguments

- input.list:

  structured input `list` with at least one element:

  - `link.functions`: `list` of functions to convert the ETAS parameters
    from the INLA scale to the ETAS scale

- magnitude:

  Magnitude of the event for which the triggering function is
  calculated, `scalar` (`default = 4`).

- n.samp:

  Number of posterior samples, `integer` (`default = 10`).

- t.end:

  Upper bound of the x-axis, `scalar` (`default = 1`).

- n.breaks:

  Number of points between 0 and `t.end` to calculate the function,
  `integer` (`default = 100`)

## Value

`ggplot` object with grey lines representing the triggering function for
each posterior sample. Black lines representing the 0.025 and 0.975
quantiles of the function values calculated for each posterior sample.
Horizontal red lines represents the 0.025 and 0.975 quantiles of the
sampled background rates.
