# Triggering function plot from posterior samples

Function to plot the ETAS triggering function corresponding to different
posterior samples

## Usage

``` r
triggering_fun_plot(
  input.list,
  post.samp = NULL,
  n.samp = 10,
  magnitude = 4,
  t.end = 1,
  n.breaks = 100
)
```

## Arguments

- input.list:

  structured input `list` with at least two elements:

  - `model.fit`: `bru` object used to sample the posterior of the ETAS
    parameters

  - `link.functions`: `list` of functions to convert the ETAS parameters
    from the INLA scale to the ETAS scale

- post.samp:

  `data.frame` containing posterior samples of the parameters. If
  `NULL`, then `n.samp` samples are generated. If `n.samp` is different
  from `nrow(post.samp)` then `n.samp` rows are uniformly sampled from
  `post.samp`. Default is `NULL`

- n.samp:

  Number of posterior samples, `integer` (`default = 10`).

- magnitude:

  Magnitude of the event for which the triggering function is
  calculated, `scalar` (`default = 4`).

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
