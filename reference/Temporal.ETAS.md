# Function to fit Hawkes process model

function to fit a temporal ETAS model using `inlabru`.

## Usage

``` r
Temporal.ETAS(
  total.data,
  M0,
  T1,
  T2,
  link.functions = NULL,
  coef.t.,
  delta.t.,
  N.max.,
  bru.opt
)
```

## Arguments

- total.data:

  Observed events: `data.frame` with columns time (ts), magnitude
  (magnitudes), event identifier (idx.p). Column names must not be
  changed.

- M0:

  Minimum magnitude threshold, `scalar`

- T1:

  Start of temporal model domain, `scalar`
  `[measure unit of sample.s$ts]`.

- T2:

  End of temporal model domain, `scalar`
  `[measure unit of sample.s$ts]`.

- link.functions:

  Functions to transform the parameters from the internal INLA scale to
  the ETAS scale. It must be a `list` of functions with names (mu, K,
  alpha, c\_, p)

- coef.t.:

  TimeBinning parameter: parameter regulating the relative length of
  successive bins, `scalar`.

- delta.t.:

  TimeBinning parameter: parameter regulating the bins' width, `scalar`.

- N.max.:

  TimeBinning parameter: parameter regulating the Number of bins (=
  `N.max` + 2), `scalar`.

- bru.opt:

  Runtime options for inlabru: See
  https://inlabru-org.github.io/inlabru/reference/bru_call_options.html,
  `list`

## Value

The fitted model as a 'bru' object, which is a list
