#' Copula transformation from a standard Normal distribution to a Gamma distribution
#'
#' @param x values from a standard Normal distribution, `vector`.
#' @param a shape parameter of the gamma distribution `scalar`.
#' @param b rate parameter of the gamma distribution `scalar`.
#'
#' @return values from a Gamma distribution with shape `a` and rate `b`, `vector` same length as `x`.
#' @export gamma.t
gamma.t <- function(x, a, b){
  bru_forward_transformation(qgamma, x, a, b)
}

#' Copula tranformation from a standard Normal distribution to a Uniform distribution
#'
#' @param x values from a standard Normal distribution, `vector`.
#' @param a minimum value for the Uniform distribution, `scalar`.
#' @param b maximum value for the Uniform distribution, `scalar`.
#'
#' @return values from a Uniform distribution between `a` and `b`, `vector` same length as `x`.
#' @export
unif.t <- function(x, a, b){
  bru_forward_transformation(qunif, x, min = a, max = b)
}

#' Copula tranformation from a standard Normal distribution to a Log-Normal distribution
#'
#' @param x values from a standard Normal distribution, `vector`.
#' @param m mean of the logarithm of the Log-Normal distribution, `scalar`.
#' @param s standard deviation of the logarithm of the Log-Normal distribution, `scalar`.
#'
#' @return Values from a Log-Normal distribution with logarithmic mean `m` and standard deviation `s`, `vector` same length as `x`.
#' @export
loggaus.t <- function(x, m, s){
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

#' Copula tranformation from a standard Normal distribution to an Exponential distribution
#'
#' @param x values from a standard Normal distribution, `vector`.
#' @param r rate of the exponential distribution, `scalar`.
#'
#' @return values from an Exponential distribution with rate `r`, `vector` same length as `x`.
#' @export exp.t
exp.t <- function(x, r){
  bru_forward_transformation(qexp, x, r)
}

#' Copula tranformation from an Exponential to a standard Normal distribution
#'
#' @param x values from an Exponential distribution, `vector`.
#' @param rate rate of the Exponential distribution, `scalar`.
#'
#' @return values from a standard Normal distribution, `vector` same length as `x`
#' @export
inv.exp.t <- function(x, rate){
  qnorm(pexp(x, rate))
}


#' Copula tranformation from an Gamma to a standard Normal distribution
#'
#' @param x values from a Gamma distribution, `vector`.
#' @param a shape parameter of the Gamma distribution, `scalar`.
#' @param b rate parameter of the Gamma distribution, `scalar`.
#'
#' @return values from a standard Normal distribution, `vector` same length as `x`
#' @export
inv.gamma.t <- function(x, a, b){
  qnorm(pgamma(x, a, b))
}

#' Copula tranformation from an Uniform to a standard Normal distribution
#'
#' @param x values from an Uniform distribution, `vector`.
#' @param a minimum of the Uniform distribution, `scalar`.
#' @param b maximum of the Uniform distribution, `scalar`.
#'
#' @return values from a standard Normal distribution, `vector` same length as `x`
#' @export
inv.unif.t <- function(x, a, b){
  qnorm(punif(x, a, b))
}

#' Copula tranformation from an Log-Normal to a standard Normal distribution
#'
#' @param x values from a Log-Normal distribution, `vector`.
#' @param m mean of the logarithm of the Log-Normal distribution, `scalar`.
#' @param s standard deviation of the logarithm of the Log-Normal distribution, `scalar`.
#'
#' @return values from a standard Normal distribution, `vector` same length as `x`
#' @export
inv.loggaus.t <- function(x, m, s){
  qnorm(plnorm(x, meanlog = m, sdlog = s))
}
