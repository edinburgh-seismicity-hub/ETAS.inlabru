# Create input list for ETAS Hawkes temporal model with catalogue

Create input list for ETAS Hawkes temporal model with catalogue

## Usage

``` r
create_input_list_temporal_withCatalogue(input_path, num.threads = NULL)
```

## Arguments

- input_path:

  path of the `txt` file containing experiment's information

- num.threads:

  Optional argument for the number of threads to be used by parallel
  processing in inlabru/INLA

## Value

The formatted input.list with the elements required for the temporal
Hawkes model
