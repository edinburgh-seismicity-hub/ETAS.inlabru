# Sample from the posterior of the ETAS parameters

Sample from the posterior of the ETAS parameters

## Usage

``` r
post_sampling(input.list, n.samp, max.batch = 1000, ncore = NULL)
```

## Arguments

- input.list:

  structured input `list` with at least two elements:

  - `model.fit`: `bru` object used to sample the posterior of the ETAS
    parameters

  - `link.functions`: `list` of functions to convert the ETAS parameters
    from the INLA scale to the ETAS scale

- n.samp:

  The number of samples to draw from the posteriors

- max.batch:

  Maximum number of posterior samples to be generated simultaneously.
  Default is 1000.

- ncore:

  Deprecated argument for controlling parallelism. Use
  `future::plan(future::multisession, workers = ncore)` (or similar) to
  configure parallelism in your code instead.

## Value

`data.frame` of posterior samples with `nrow = n.samp` and columns
`mu, K, alpha, c, p` corresponding to ETAS parameters.
