
###################################
## functions to run automatically
## Conversion of ETAS parameters to internal scale and back
## MN:: Check this is the only use...

### Forward link function

#' Gamma copula transformation: Conversion of ETAS para to internal scale
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

#
#' Uniform copula transformation: Conversion of ETAS para to internal scale
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

#' Log-gaussian copula transformation: Conversion of ETAS para to internal scale
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


#' Inverse exponential link function:
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


#' integrated triggering function - used by Inlabru ALSO USED IN GENERATION SAMPLES
#'
#' @param theta ETAS parameters `data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`.
#' @param ti Time of parent event.
#' @param mi Magnitude of parent event
#' @param M0 Minimum magnitude threshold
#' @param T1 Start of temporal model domain.
#' @param T2 End of temporal model domain.
#'
#' @return
#' @export
#'
#' @examples
log.Lambda_h2 <- function(theta, ti, mi, M0, T1, T2){
  th <- theta
  T.low <- pmax(T1, ti)#sapply(ti, \(x) max(T1, x))

  gamma.l <- (T.low - ti)/th$c
  gamma.u <- (T2 - ti)/th$c
  w.l <- (1 + gamma.l)^(1-th$p)
  w.u <- (1 + gamma.u)^(1-th$p)
  # output
  log(th$K) + th$alpha*(mi - M0) + log(th$c) - log(th$p - 1) + log1p(w.l - 1) + log1p(-w.u/w.l)
}

#' Integrated ETAS time-triggering function
#'
#' @param theta ETAS parameters data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)
#' @param th Time of past event? [days]
#' @param T2
#'
#' @return
#' @export
#'
#' @examples
# It <- function(theta, th, T2){
Int.ETAS.time.trig.function <- function(theta, th, T2){
  gamma.u <- (T2 - th)/theta$c
  ( theta$c/(theta$p - 1) )*(1 - (gamma.u + 1)^( 1-theta$p) )
}


#' Inverse of integrated ETAS time-triggering function
#'
#' @param theta ETAS parameters data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)
#' @param omega
#' @param th
#'
#' @return
#' @export
#'
#' @examples
#Inv.It <- function(theta, omega, th){
Inv.Int.ETAS.time.trig.function <- function(theta, omega, th){
  th + theta$cc*( ( 1 - omega * (1/theta$c)*(theta$p - 1) )^( -1/(theta$p - 1) ) - 1)
}
