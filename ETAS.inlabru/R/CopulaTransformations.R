
# gamma copula transformation
#' Title
#'
#' @param x
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
gamma.t <- function(x, a, b){
  bru_forward_transformation(qgamma, x, a, b)
}
# uniform copula transformation
#' Title
#'
#' @param x
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
unif.t <- function(x, a, b){
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
#' Title
#'
#' @param x
#' @param m
#' @param s
#'
#' @return
#' @export
#'
#' @examples
loggaus.t <- function(x, m, s){
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

#' Title
#'
#' @param x
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
inv.exp.t <- function(x, rate){
  qnorm(pexp(x, rate))
}


#' Inverse gamma copula transformation:
#'
#' @param x
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
inv.gamma.t <- function(x, a, b){
  qnorm(pgamma(x, a, b))
}

#' Inverse uniform copula transformation:
#'
#' @param x
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
inv.unif.t <- function(x, a, b){
  qnorm(punif(x, a, b))
}

#' Inverse log-gaussian copula transformation:
#'
#' @param x
#' @param m
#' @param s
#'
#' @return
#' @export
#'
#' @examples
inv.loggaus.t <- function(x, m, s){
  qnorm(plnorm(x, meanlog = m, sdlog = s))
}
