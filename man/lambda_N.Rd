% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysingPosteriors.R
\name{lambda_N}
\alias{lambda_N}
\title{Calculate the integral of the ETAS conditional intensity}
\usage{
lambda_N(th.mu, th.K, th.alpha, th.c, th.p, T1, T2, M0, Ht, link.functions)
}
\arguments{
\item{th.mu}{Background rate, \code{mu}, on the internal parameter scale}

\item{th.K}{ETAS triggering parameter \code{K} on the internal parameter scale}

\item{th.alpha}{ETAS triggering parameter \code{alpha} on the internal parameter scale}

\item{th.c}{ETAS triggering parameter \code{c} on the internal parameter scale}

\item{th.p}{ETAS triggering parameter \code{p} on the internal parameter scale}

\item{T1}{Start of temporal model domain.}

\item{T2}{End of temporal model domain.}

\item{M0}{Minimum magnitude threshold}

\item{Ht}{History of the process, or set of known events in the interval. It must be a \code{data.frame} with columns \code{ts} (time) and \code{magnitudes} (magnitudes).}

\item{link.functions}{\code{list} of functions to transform the parameters from the internal scale to the ETAS scale}
}
\value{
Integral of the ETAS conditional intensity between \code{T1} and \code{T2} with minimum magnitude \code{M0},
i.e. the expected number of events.
}
\description{
Calculate the number of events in a time interval T1 to T2 given imposed events and ETAS
}
