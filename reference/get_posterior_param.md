# Retrieve posterior distribution of ETAS parameters

Retrieve posterior distribution of ETAS parameters

## Usage

``` r
get_posterior_param(input.list)
```

## Arguments

- input.list:

  input.list structured input `list` with at least two elements:

  - `model.fit`: `bru` object used to sample the posterior of the ETAS
    parameters

  - `link.functions`: `list` of functions to convert the ETAS parameters
    from the INLA scale to the ETAS scale

## Value

A `list` of two elements:

- `post.df` : `data.frame` with the posterior distributions of the
  parameters with columns `x` (value of the parameter), `y` (value of
  the posterior), `param` (parameter name)

- `post.plot` : `ggplot` object showing the posterior distribution of
  each parameter
