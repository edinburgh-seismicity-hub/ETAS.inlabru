# 2b Temporal Model: Presenting samples drawn from the ETAS priors

``` r
library(ETAS.inlabru)
library(ggplot2)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
# library(tidyquant)
```

## Sampling ETAS Priors

It is important to check that the priors are broad enough that we expect
the the model parameterisation to lie within them. This notebook shows
how to draw samples from the priors and present them and the resulting
triggering functions.
