% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateSyntheticCatalogues.R
\name{sample_GR_magnitudes}
\alias{sample_GR_magnitudes}
\title{Return a sample of magnitudes drawn from the GR distribution}
\usage{
sample_GR_magnitudes(n, beta.p, M0)
}
\arguments{
\item{n}{Number of events in the sample.}

\item{beta.p}{Related to the b-value via \verb{b ln(10)}.}

\item{M0}{Minimum magnitude for the sample.}
}
\value{
A list of magnitudes of length \code{n} drawn from a GR distribution.
}
\description{
Return a sample of magnitudes drawn from the GR distribution
}
\examples{
sample_GR_magnitudes(n = 100, beta.p = log(10), M0 = 2.5)

}
