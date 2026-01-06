# ETAS.inlabru: Bayesian ETAS model for modelling seismic sequences with inlabru

Modelling and inversion of ETAS model of seismicity using inlabru The
Epidemic Type Aftershock Sequence (ETAS) model is designed to model
earthquakes that are triggered by previous events. In statistics, this
is referred to as a Hawkes process. The code can be used to generate
synthetic ETAS catalogues which can also include some seeded events to
model specific sequences. We also implement a Bayesian inversion scheme
using the Integrated Nested Laplace Approximation (INLA) using inlabru.
For the temporal model, given a training catalogue of times and
magnitudes, the code returns the joint posteriors for all the ETAS
parameters. In the future roadmap, we will include tools to model the
spatial distribution and spatio-temporal evolution of seismic sequences.

## See also

Useful links:

- <https://edinburgh-seismicity-hub.github.io/ETAS.inlabru/>

- <https://github.com/edinburgh-seismicity-hub/ETAS.inlabru>

- Report bugs at
  <https://github.com/edinburgh-seismicity-hub/ETAS.inlabru/issues>

## Author

**Maintainer**: Francesco Serafini <francesco.serafini@newcastle.ac.uk>
([ORCID](https://orcid.org/0000-0003-0154-6200))

Authors:

- Mark Naylor <mark.naylor@ed.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-3761-5522)) \[thesis advisor\]

- Finn Lindgren <Finn.Lindgren@ed.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-5833-2011)) \[thesis advisor\]
