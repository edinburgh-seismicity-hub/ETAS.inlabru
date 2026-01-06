# Inverse of integrated Omori's law

Inverse of integrated Omori's law

## Usage

``` r
Inv_Int_ETAS_time_trig_function(theta, omega, th)
```

## Arguments

- theta:

  ETAS parameters `list(mu=mu, K=K, alpha=alpha, c=c, p=p)`

- omega:

  Value of the integral to be inverted, `vector`

- th:

  Time from which the integral is calculated `scalar`

## Value

Value of the start of the temporal domain used to calculate the integral

## Details

Considering the integral of Omori's law \$\$\omega =
\int\_{t_h}^{T_2}\left(\frac{t - t_h}{c} + 1\right)^{-p} dt\$\$ The
function applied to the value \\\omega\\ returns the value of \\t_h\\.
