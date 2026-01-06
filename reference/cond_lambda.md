# ETAS conditional intensity - used by `inlabru`

Function to calculate the value of the ETAS model conditional intensity
at a specified time given the history of the process.

## Usage

``` r
cond_lambda(theta, t, th, mh, M0)
```

## Arguments

- theta:

  ETAS parameters as `list` with names `mu`, `K`, `alpha`, `c`, `p`

- t:

  Time at which the conditional intensity is evaluated, `scalar`

- th:

  Time of the events in the history of the process, `vector`

- mh:

  Magnitudes of the events in the history of the process, `vector`

- M0:

  Minimum magnitude threshold

## Value

Value of the ETAS conditional intensity calculated at `t` with history
`th, mh`, `scalar`

## Details

The function takes a single value `t` and returns the ETAS conditional
intensity calculated at `t` with history `th, mh`. The ETAS conditional
intensity is given by \$\$\lambda(t \| \mathcal H_t) = \mu + \sum\_{h:
(t_h,m_h) \in \mathcal H_t} K e^{\alpha(m_h - M_0)} \left( \frac{t -
t_h}{c} + 1\right)^{-p}\$\$

Do not use when `t` is a vector.
