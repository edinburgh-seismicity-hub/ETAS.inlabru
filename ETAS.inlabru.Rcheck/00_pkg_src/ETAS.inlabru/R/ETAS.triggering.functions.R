
#' integrated triggering function - used by Inlabru ALSO USED IN GENERATION SAMPLES
#'
#' @param theta data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)
#' @param ti
#' @param mi
#' @param M0
#' @param T1
#' @param T2
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
#' @param theta
#' @param th
#' @param T2
#'
#' @return
#' @export
#'
#' @examples
# It <- function(theta, th, T2){
Int.ETAS.time.trig.function <- function(theta, th, T2){
  gamma.u <- (T2 - th)/theta$c
  ( theta$c/(theta$p - 1) )*(1 - (gamma.u + 1)^(1-theta$p) )
}


#' Inverse of integrated ETAS time-triggering function
#'
#' @param theta
#' @param omega
#' @param th
#'
#' @return
#' @export
#'
#' @examples
#Inv.It <- function(theta, omega, th){
Inv.Int.ETAS.time.trig.function <- function(theta, omega, th){
  th + theta$c*( ( 1 - omega * (1/theta$c)*(theta$p - 1) )^( -1/(theta$p - 1) ) - 1)
}
