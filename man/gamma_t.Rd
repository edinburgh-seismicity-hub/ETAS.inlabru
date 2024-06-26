% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CopulaTransformations.R
\name{gamma_t}
\alias{gamma_t}
\title{Copula transformation from a standard Normal distribution to a Gamma distribution}
\usage{
gamma_t(x, a, b)
}
\arguments{
\item{x}{values from a standard Normal distribution, \code{vector}.}

\item{a}{shape parameter of the gamma distribution \code{scalar}.}

\item{b}{rate parameter of the gamma distribution \code{scalar}.}
}
\value{
values from a Gamma distribution with shape \code{a} and rate \code{b}, \code{vector} same length as \code{x}.
}
\description{
Copula transformation from a standard Normal distribution to a Gamma distribution
}
\seealso{
Other copula-transformations: 
\code{\link{exp_t}()},
\code{\link{inv_exp_t}()},
\code{\link{inv_gamma_t}()},
\code{\link{inv_loggaus_t}()},
\code{\link{inv_unif_t}()},
\code{\link{loggaus_t}()},
\code{\link{unif_t}()}
}
\concept{copula-transformations}
