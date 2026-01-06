# Title

Title

## Usage

``` r
sample_temporal_injection_events(
  a = 50,
  V.i = 1,
  tau = 10,
  beta.p,
  M0,
  T.i,
  T2
)
```

## Arguments

- a:

  Induced event rate per unit volume.

- V.i:

  Injected volume

- tau:

  Decay rate `[days]`.

- beta.p:

  Related to the b-value via `b ln(10)`.

- M0:

  Minimum magnitude threshold.

- T.i:

  Time of injection `[days]`.

- T2:

  End of temporal model domain `[days]`.

## Value

Catalogue of parent events induced by injection;
`data.frame(times, magnitudes)`
