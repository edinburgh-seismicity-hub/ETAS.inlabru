# Injection Rate function calculations

Forward time integrated function for exponential rate decay, and its
inverse

## Usage

``` r
IntInjectionIntensity(a = 50, V.i = 1, tau = 10, T.i, T2)

Inv_IntInjectionIntensity(
  a = 50,
  V.i = 1,
  tau = 10,
  T.i,
  number.injected.events
)
```

## Arguments

- a:

  Event rate per unit volume injected

- V.i:

  Injected volume

- tau:

  Decay rate `[days]`

- T.i:

  Time of injection event

- T2:

  End of temporal model domain

- number.injected.events:

  The number of expected injected events, used for the inverse.

## Value

`IntInjectionIntensity` returns the forward time integrated function for
exponential rate decay.

`Inv_IntInjectionIntensity` returns the end time corresponding to a
given expected number of injected events.
