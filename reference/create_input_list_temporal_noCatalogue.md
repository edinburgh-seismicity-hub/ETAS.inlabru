# Create input list for ETAS Hawkes temporal model without catalogue

Create input list for ETAS Hawkes temporal model without catalogue

## Usage

``` r
create_input_list_temporal_noCatalogue(input_path, num.threads = NULL)
```

## Arguments

- input_path:

  Input file and path as a string

- num.threads:

  Optional argument for the number of threads to be used by parallel
  processing by inlabru/INLA

## Value

The formatted input.list with the elements required for the temporal
Hawkes model

## Examples

``` r
create_input_list_temporal_noCatalogue(
  system.file(
    "extdata",
    "user_input_synthetic_noCatalogue.txt",
    package = "ETAS.inlabru"
  )
)
#> $catalog
#> NULL
#> 
#> $catalog.bru
#> NULL
#> 
#> $time.int
#> NULL
#> 
#> $T12
#> [1] "T1"  " T2"
#> 
#> $lat.int
#> [1] -90  90
#> 
#> $lon.int
#> [1] -180  180
#> 
#> $M0
#> NULL
#> 
#> $mu.init
#> [1] 0.25
#> 
#> $K.init
#> [1] 0.3
#> 
#> $alpha.init
#> [1] 1.6
#> 
#> $c.init
#> [1] 0.2
#> 
#> $p.init
#> [1] 1.1
#> 
#> $a_mu
#> [1] 0.5
#> 
#> $b_mu
#> [1] 0.5
#> 
#> $a_K
#> [1] -1
#> 
#> $b_K
#> [1] 0.5
#> 
#> $a_alpha
#> [1] 0
#> 
#> $b_alpha
#> [1] 10
#> 
#> $a_c
#> [1] 0
#> 
#> $b_c
#> [1] 1
#> 
#> $a_p
#> [1] 1
#> 
#> $b_p
#> [1] 2
#> 
#> $max_iter
#> [1] 100
#> 
#> $max_step
#> NULL
#> 
#> $link.functions
#> $link.functions$mu
#> function (x) 
#> gamma_t(x, a_mu, b_mu)
#> <bytecode: 0x56031f501570>
#> <environment: 0x56031f4fed58>
#> 
#> $link.functions$K
#> function (x) 
#> loggaus_t(x, a_K, b_K)
#> <bytecode: 0x56031f4fda18>
#> <environment: 0x56031f4fed58>
#> 
#> $link.functions$alpha
#> function (x) 
#> unif_t(x, a_alpha, b_alpha)
#> <bytecode: 0x56031f4fdcf0>
#> <environment: 0x56031f4fed58>
#> 
#> $link.functions$c_
#> function (x) 
#> unif_t(x, a_c, b_c)
#> <bytecode: 0x56031f4fe000>
#> <environment: 0x56031f4fed58>
#> 
#> $link.functions$p
#> function (x) 
#> unif_t(x, a_p, b_p)
#> <bytecode: 0x56031f4fe2d8>
#> <environment: 0x56031f4fed58>
#> 
#> 
#> $bru.opt.list
#> $bru.opt.list$bru_verbose
#> [1] 3
#> 
#> $bru.opt.list$bru_max_iter
#> [1] 100
#> 
#> $bru.opt.list$num.threads
#> NULL
#> 
#> $bru.opt.list$bru_initial
#> $bru.opt.list$bru_initial$th.mu
#> [1] -0.2978078
#> 
#> $bru.opt.list$bru_initial$th.K
#> [1] -0.4079456
#> 
#> $bru.opt.list$bru_initial$th.alpha
#> [1] -0.9944579
#> 
#> $bru.opt.list$bru_initial$th.c
#> [1] -0.8416212
#> 
#> $bru.opt.list$bru_initial$th.p
#> [1] -1.281552
#> 
#> 
#> 
#> $coef.t
#> [1] 1
#> 
#> $delta.t
#> [1] 0.1
#> 
#> $Nmax
#> [1] 8
#> 
#> $n.periods
#> [1] 120
#> 
#> $period.length
#> [1] 1
#> 
#> $start.date.fore
#> NULL
#> 
#> $magnitude.update
#> [1] 5.5
#> 
#> $output.name
#> [1] "report_ETAS"
#> 
```
