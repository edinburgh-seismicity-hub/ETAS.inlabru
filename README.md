# ETAS.inlabru
R package that implements the ETAS Hawkes process for modelling seismicity

# Authors

- Francesco Serafini
- Dr Mark Naylor, School of GeoSciences, University of Edinburgh
- Prof Finn Lindgren, School of Mathematics, University of Edinburgh
- Dr Kirsty Bayliss, Global Earthquake Model (GEM)

# Installation

For ETAS.inlabru to work, we need to install both R-INLA and inlabru:

For inlabru (see https://github.com/inlabru-org/inlabru) :
`remotes::install_github("inlabru-org/inlabru", ref="devel")`

For R-INLA (see https://www.r-inla.org/download-install) :
`install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)`


# Terminology and planning suggestions

## File structure in package

### Data files

### Code files

- ETAS.triggering.function.R : Contains the ETAS specific model functions
- HawkesProcess.R : Generic Hawkes code that is intended for integration back into inlabru
- generateSyntheticCatalogues.R : Contains the iterative Hawkes functions for generating triggered events but the actual triggering functions reside in the ETAS file above so we could introduce other models here
- temporalBinning.R : Code to generate time bins to make integration scheme efficient
- plottingFunctions.R : lets put all the standard plotting functions into here
- setupInlabruInputs.R : Put the functions for generating input.list in here

## Terminology

- Let's be specific about when we are doing just temporal so we have clear function names for the spatial and spatial-temporal later on
- I changed the theta from c(mu, K, ...) to a  df <- data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p). This means we should refer to the values as df$mu which is unambiguous etc
  - Might this do anything bad?
- Sometimes you used th for theta and sometimes for historic times...
  - I have tried to modify to just theta

# What has been done

## Implemented

- Generation of synthetic ETAS catalogues by ETAS.inlabru with a demonstration in the notebook

## In development

- Add inversion modelling based on original code
- Modify the implementation so that the generic Hawkes code can go into inlabru and the ETAS triggering function code stay in this package

## Roadmap

- Implement incompleteness fix
- Implement pre-model domain history
- Integrate spatial modelling
