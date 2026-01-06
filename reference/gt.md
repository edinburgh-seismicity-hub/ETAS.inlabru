# ETAS triggering function - used by `inlabru`

The function returns the value of the ETAS triggering function at a
specified time `t` for the points in the history `th, mh`

## Usage

``` r
gt(theta, t, th, mh, M0)
```

## Arguments

- theta:

  ETAS parameters as `list` with names `K`, `alpha`, `c`, `p`

- t:

  Time at which the function is calculated, `scalar` or `vector`

- th:

  Time of events in the history in `[days, months,...]`, `scalar` or
  `vector`

- mh:

  Magnitude of events in the history, `scalar` or `vector`

- M0:

  Minimum magnitude threshold, `scalar`

## Value

value of the ETAS triggering function evaluated at `t` with history
`th`, `mh`.

## Details

The ETAS triggering function to be evaluated is \$\$g(t - t_h \| m_h) =
K e^{\alpha(m_h - M_0)} \left( \frac{t - t_h}{c} + 1\right)^{-p}\$\$
Where \\K, \alpha, c \> 0\\, and \\p \geq 1\\ are the ETAS parameters,
\\t\\ (`t`) is the time at which the function is evaluated, considering
the past observation \\t_h, m_h\\ (`th, mh`). The function returns 0 if
\\t_h \> t\\. If \\t\\ is a scalar and \\t_h, m_h\\ are vectors than the
function returns a vector, same if \\t\\ is a vector and \\t_h, m_h\\
are scalars, or if \\t, t_h, m_h\\ are vectors of the same length.

Do not use if \\t\\ and \\t_h, m_h\\ are vectors of different lengths.
