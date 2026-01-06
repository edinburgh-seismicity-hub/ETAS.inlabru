# Function to calculate Omori's law

Function to calculate Omori's law

## Usage

``` r
omori(theta, t, ti)
```

## Arguments

- theta:

  ETAS parameters (`list(mu = mu, K = K, alpha = alpha, c = c, p = p`),
  only parameters `c` and `p` are used

- t:

  Time at which Omori's law is evaluated

- ti:

  Time of the event in the history

## Value

Value of Omori's law at point `t` for and event happened in `ti`
