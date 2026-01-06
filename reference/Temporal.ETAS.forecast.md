# Title

Title

## Usage

``` r
Temporal.ETAS.forecast(post.samp, n.cat, beta.p, M0, T1, T2, Ht, ncore = NULL)
```

## Arguments

- post.samp:

  a `data.frame` containing samples from the posterior distribution of
  ETAS parameters. Each row of the `data.frame` corresponds to a
  different sample and the parameters are in the order \\\mu\\, \\K\\,
  \\\alpha\\, \\c\\, \\p\\.

- n.cat:

  number of synthetic catalogues composing the forecast. If `n.cat` is
  greater than `nrow(post.samp)`, then, `n.cat` rows are sampled
  uniformly and with replacement from the `post.samp`. If `n.cat` is
  smaller than `nrow(post.samp)`, then, `n.cat` rows are sampled
  uniformly and without replacement from the `post.samp`. If `n.cat` is
  `NULL` or equal to `nrow(post.samp)`, `post.samp` is used at it is and
  `nrow(post.samp)` catalogues are generated.

- beta.p:

  parameter of the magnitude distribution

- M0:

  cutoff magnitude, all the synthetic events will have magnitude greater
  than this value.

- T1:

  starting time of the forecast

- T2:

  end time of the forecast

- Ht:

  set of known events

- ncore:

  Deprecated argument for controlling parallelism. Use
  `future::plan(future::multisession, workers = ncore)` (or similar) to
  configure parallelism in your code instead.

## Value

a `list` with two elements: `fore.df` is a `data.frame` containing all
the synthetic catalogues composing the forecast. The `data.frame` has
four columns, `ts` for the occurrence time, `magnitudes` for the
magnitude, `gen` with the generation of the event, and `cat.idx` with
the catalogue identifier The second element of the output `list` is
`n.cat` which is the number of synthetic catalogues generated.

## See also

[`generate_temporal_ETAS_synthetic()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/generate_temporal_ETAS_synthetic.md)
