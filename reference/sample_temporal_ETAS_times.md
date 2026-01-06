# Sampling times for events triggered by a parent at th according to the ETAS triggering function

Sampling times for events triggered by a parent at th according to the
ETAS triggering function

## Usage

``` r
sample_temporal_ETAS_times(theta, n.ev, th, T2)
```

## Arguments

- theta:

  ETAS parameters `list(mu=mu, K=K, alpha=alpha, c=c, p=p)`.

- n.ev:

  Number of events to return in the sample in time domain `(th, T2]`.

- th:

  Time of the parent event producing `n.ev` daughters.

- T2:

  End time of model domain.

## Value

t.sample A list of times in the interval `[0, T2]` distributed according
to the ETAS triggering function.
