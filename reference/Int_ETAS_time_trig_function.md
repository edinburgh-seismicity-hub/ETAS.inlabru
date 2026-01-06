# Integrated Omori's law

Integrated Omori's law

## Usage

``` r
Int_ETAS_time_trig_function(theta, th, T2)
```

## Arguments

- theta:

  ETAS parameters `list(mu=mu, K=K, alpha=alpha, c=c, p=p)`

- th:

  Time of past event `[days]` and start of temporal domain, `vector`.

- T2:

  End of temporal domain, `scalar`.

## Value

Value of the integral of Omori's law

## Details

The function returns the integral of Omori's law, namely
\$\$\int\_{t_h}^{T_2} \left(\frac{t - t_h}{c} + 1\right)^{-p} dt\$\$
