# Function to plot Omori's law corresponding to different posterior samples

Function to plot Omori's law corresponding to different posterior
samples

## Usage

``` r
omori_plot_posterior(
  input.list,
  post.samp = NULL,
  n.samp = 10,
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
  `post.samp`. Default is `NULL`.

- n.samp:

  Number of posterior samples, `integer` (`default = 10`).

- t.end:

  Upper bound of the x-axis, `scalar` (`default = 1`).

- n.breaks:

  Number of points between 0 and `t.end` to calculate the function,
  `integer` (`default = 100`).

## Value

A ggplot object

## See also

[`create_input_list_temporal_noCatalogue()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/create_input_list_temporal_noCatalogue.md),
[`create_input_list_temporal_withCatalogue()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/create_input_list_temporal_withCatalogue.md)
