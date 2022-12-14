% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CopulaTransformations.R
\name{inv.loggaus.t}
\alias{inv.loggaus.t}
\title{Copula tranformation from an Log-Normal to a standard Normal distribution}
\usage{
inv.loggaus.t(x, m, s)
}
\arguments{
\item{x}{values from a Log-Normal distribution, \code{vector}.}

\item{m}{mean of the logarithm of the Log-Normal distribution, \code{scalar}.}

\item{s}{standard deviation of the logarithm of the Log-Normal distribution, \code{scalar}.}
}
\value{
values from a standard Normal distribution, \code{vector} same length as \code{x}
}
\description{
Copula tranformation from an Log-Normal to a standard Normal distribution
}
