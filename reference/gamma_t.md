# Copula transformation from a standard Normal distribution to a Gamma distribution

Copula transformation from a standard Normal distribution to a Gamma
distribution

## Usage

``` r
gamma_t(x, a, b)
```

## Arguments

- x:

  values from a standard Normal distribution, `vector`.

- a:

  shape parameter of the gamma distribution `scalar`.

- b:

  rate parameter of the gamma distribution `scalar`.

## Value

values from a Gamma distribution with shape `a` and rate `b`, `vector`
same length as `x`.

## See also

Other copula-transformations:
[`exp_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/exp_t.md),
[`inv_exp_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/inv_exp_t.md),
[`inv_gamma_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/inv_gamma_t.md),
[`inv_loggaus_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/inv_loggaus_t.md),
[`inv_unif_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/inv_unif_t.md),
[`loggaus_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/loggaus_t.md),
[`unif_t()`](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/unif_t.md)
