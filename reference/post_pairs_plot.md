# Plot the posterior densities of the ETAS parameters

Plot the posterior densities of the ETAS parameters

## Usage

``` r
post_pairs_plot(
  input.list = NULL,
  n.samp = NULL,
  post.samp = NULL,
  max.batch = 1000
)
```

## Arguments

- input.list:

  structured input `list` with at least two elements:

  - `model.fit`: `bru` object used to sample the posterior of the ETAS
    parameters

  - `link.functions`: `list` of functions to convert the ETAS parameters
    from the INLA scale to the ETAS scale

- n.samp:

  The number of samples to draw from the posteriors for the plot

- post.samp:

  `data.frame` with columns mu, K, alpha, c, p and rows corresponding to
  different posterior samples. When `NULL` the function samples the
  joint posterior distribution `n.samp` times. The default is `NULL`.

- max.batch:

  parameter of
  [post_sampling](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/post_sampling.md)
  function to be used in case `post.samp = NULL`

## Value

`list`: with elements

- `post.samp.df`:`data.frame` of posterior samples with `nrow = n.samp`
  and columns `mu, K, alpha, c, p` corresponding to ETAS parameters. If
  `post.samp` is not `NULL` it returns `post.samp`

- `pair.plot`: `ggplot` object reporting the pair plot between
  parameters samples. It is obtained using the `ggpairs` function of the
  `Ggally` library
