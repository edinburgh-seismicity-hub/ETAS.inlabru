# Fit temporal ETAS model

Fits the remporal ETAS model and returns the results. This function
decomposes the input.list for the \`Hawkes.bru2“ function.

## Usage

``` r
Temporal.ETAS.fit(input.list)
```

## Arguments

- input.list:

  All input data and parameters are passed to inlabru via this
  structured `list`. This is the output of the function
  [create_input_list_temporal_withCatalogue](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/create_input_list_temporal_withCatalogue.md)
  or
  [create_input_list_temporal_noCatalogue](https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/reference/create_input_list_temporal_noCatalogue.md)

## Value

The fitted model as a `bru` object, which is a list
