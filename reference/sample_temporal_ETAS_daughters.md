# Sample daughter events from one parent using the ETAS model

Generate a sample of new events `data.frame(t_i, M_i)` of length `n.ev`
for one parent event occuring at time `t_h` using the ETAS model.

## Usage

``` r
sample_temporal_ETAS_daughters(theta, beta.p, th, n.ev, M0, T1, T2)
```

## Arguments

- theta:

  ETAS parameters `list(mu=mu, K=K, alpha=alpha, c=c, p=p)`.

- beta.p:

  Slope of GR relation: beta = b ln(10).

- th:

  Time of parent event `[days]`.

- n.ev:

  The number of events to be placed.

- M0:

  Minimum magnitude in synthetic catalogue.

- T1:

  Start time for synthetic catalogue `[days]`.

- T2:

  End time for synthetic catalogue `[days]`.

## Value

Generate a sample of new events `data.frame(t_i, M_i)` from one parent
