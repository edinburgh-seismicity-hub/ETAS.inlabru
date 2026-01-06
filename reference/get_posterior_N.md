# Plot the posterior distribution of the expected number of events

Plot the posterior distribution of the expected number of events

## Usage

``` r
get_posterior_N(input.list, domain.extension = 0.1)
```

## Arguments

- input.list:

  Which has combined the input file (for link functions) and bru output
  (for marginals)

- domain.extension:

  Percentage of posterior quantiles to extend the domain specified as
  `scalar`. Default is set to 0.10.

## Value

A `list` of three objects:

- `post.df`: `data.frame` containing posterior informations on the
  posterior distribution of the number of events

- `post.plot` : `ggplot` object showing the posterior distribution of
  the expected number of events

- `post.plot.shaded` : `ggplot` object showing the posterior
  distribution of the expected number of events, shaded region
  corresponds to the 0.025 and 0.975 quantiles of the distribution of
  the distribution of the number of events.
