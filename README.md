# ETAS.inlabru
R package that implements the ETAS Hawkes process for modelling seismicity

# Terminology and planning suggestions

File structure in package

- ETAS.triggering.function.R : Contains the ETAS specific model functions
- generateSyntheticCatalogues.R : Contains the iterative Hawkes functions for generating triggered events but the actual triggering functions reside in the ETAS file above so we could introduce other models here
- plottingFunctions.R : lets put all the standard plotting functions into here
- setupInlabruInputs.R : Put the functions for generating input.list in here

Terminology

- Let's be specific about when we are doing just temporal so we have clear function names for the spatial and spatial-temporal later on
- I changed the theta from c(mu, K, ...) to a  df <- data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p). This means we should refer to the values as df$mu which is unambiguous etc
  - Might this do anything bad?
- Sometimes you used th for theta and sometimes for historic times...
  - I have tried to modify to just theta

# What has been done

## Implemented

- Generation of synthetic ETAS catalogues by ETAS.inlabru with a demonstration in the notebook

## In development

- We need to tidy up the generation stuff

## Roadmap

