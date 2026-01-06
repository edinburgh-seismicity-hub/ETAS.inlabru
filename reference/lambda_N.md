# Calculate the integral of the ETAS conditional intensity

Calculate the number of events in a time interval T1 to T2 given imposed
events and ETAS

## Usage

``` r
lambda_N(th.mu, th.K, th.alpha, th.c, th.p, T1, T2, M0, Ht, link.functions)
```

## Arguments

- th.mu:

  Background rate, `mu`, on the internal parameter scale

- th.K:

  ETAS triggering parameter `K` on the internal parameter scale

- th.alpha:

  ETAS triggering parameter `alpha` on the internal parameter scale

- th.c:

  ETAS triggering parameter `c` on the internal parameter scale

- th.p:

  ETAS triggering parameter `p` on the internal parameter scale

- T1:

  Start of temporal model domain.

- T2:

  End of temporal model domain.

- M0:

  Minimum magnitude threshold

- Ht:

  History of the process, or set of known events in the interval. It
  must be a `data.frame` with columns `ts` (time) and `magnitudes`
  (magnitudes).

- link.functions:

  `list` of functions to transform the parameters from the internal
  scale to the ETAS scale

## Value

Integral of the ETAS conditional intensity between `T1` and `T2` with
minimum magnitude `M0`, i.e. the expected number of events.
