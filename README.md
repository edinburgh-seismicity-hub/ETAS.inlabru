
<!-- README.md is generated from README.Rmd. Please edit the .Rmd file -->

# ETAS.inlabru

<!-- badges: start -->
<!-- badges: end -->

R package that implements the ETAS Hawkes process for modelling
seismicity

Online documentation:
<https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/>

## Authors

- [Dr Francesco
  Serafini](https://scholar.google.com/citations?user=NVDOxTcAAAAJ&hl=en)
- [Dr Mark Naylor](https://blogs.ed.ac.uk/mnaylor/) , School of
  GeoSciences, University of Edinburgh
- [Prof Finn Lindgren](https://www.maths.ed.ac.uk/~flindgre/) , School
  of Mathematics, University of Edinburgh
- [Dr Kirsty
  Bayliss](https://www.linkedin.com/in/kirsty-bayliss-9a6604a1/?originalSubdomain=uk)
  , Global Earthquake Model (GEM)

## Funding

- This study was funded by yhe Real-Time Earthquake Risk Reduction for a
  Resilient Europe [RISE project](http://www.rise-eu.org/home/) , which
  has received funding from the European Union’s Horizon 2020 Research
  and Innovation Program under grant Agreement 821115.
- Naylor was additionally funded by the NSFGEO-NERC grant NE/R000794/1.
- Bayliss was funded by an EPSRC Studentship.

## Installation

For ETAS.inlabru to work, we need to install both R-INLA and inlabru:

For inlabru (see <https://inlabru-org.github.io/inlabru/>):

- CRAN release,

``` r
install.packages("inlabru")
```

- or development version,

``` r
# install.packages("remotes")
remotes::install_github("inlabru-org/inlabru")
```

For R-INLA (see <https://www.r-inla.org/download-install>):

``` r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/testing"),
  dep = TRUE
)
```

You can install the development version of ETAS.inlabru from
[GitHub](https://github.com/) with

``` r
# install.packages("remotes")
remotes::install_github("edinburgh-seismicity-hub/ETAS.inlabru")
```

## Terminology and planning suggestions

### File structure in package

- ETAS.triggering.function.R : Contains the ETAS specific model
  functions
- HawkesProcess.R : Generic Hawkes code that is intended for integration
  back into inlabru
- generateSyntheticCatalogues.R : Contains the iterative Hawkes
  functions for generating triggered events but the actual triggering
  functions reside in the ETAS file above so we could introduce other
  models here
- temporalBinning.R : Code to generate time bins to make integration
  scheme efficient
- plottingFunctions.R : lets put all the standard plotting functions
  into here
- setupInlabruInputs.R : Put the functions for generating input.list in
  here

## Terminology

- Let’s be specific about when we are doing just temporal so we have
  clear function names for the spatial and spatial-temporal later on
- I changed the theta from `c(mu, K, ...)` to a
  `df <- data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`. This means we
  should refer to the values as `df$mu` which is unambiguous etc
  - Might this do anything bad?
- Sometimes you used `th` for `theta` and sometimes for historic times…
  - I have tried to modify to just `theta`

# What has been done

## Implemented

- Generation of synthetic ETAS catalogues by `ETAS.inlabru` with a
  demonstration in the notebook

## In development

- Add inversion modelling based on original code
- Modify the implementation so that the generic Hawkes code can go into
  inlabru and the ETAS triggering function code stay in this package

## Roadmap

- Integrate spatial modelling
