# Plot Omori's law for prior samples

Function to plot Omori's law corresponding to different prior samples

## Usage

``` r
omori_plot_prior(input.list, n.samp = 10, t.end = 1, n.breaks = 100)
```

## Arguments

- input.list:

  structured input `list` with at least one element:

  - `link.functions`: `list` of functions to convert the ETAS parameters
    from the INLA scale to the ETAS scale

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
