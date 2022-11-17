

##################################################################################
## Functions for ETAS conditional intensity and integrated conditional intensity
#################################################################################


#' ETAS triggering function - used by `inlabru`
#' @description The function returns the value of the ETAS triggering function at a specified time `t` for the points in the history `th, mh`
#' @param theta ETAS parameters as `list` with names `K`, `alpha`, `c`, `p`
#' @param t Time at which the function is calculated, `scalar` or `vector`
#' @param th  Time of events in the history in [days, months,...], `scalar` or `vector`
#' @param mh Magnitude of events in the history, `scalar` or `vector`
#' @param M0 Minimum magnitude threshold, `scalar`
#'
#' @return value of the ETAS triggering function evaluated at `t` with history `th`, `mh`.
#' @export
#' @details The ETAS triggering function to be evaluated is
#' \deqn{g(t - t_h | m_h) = K e^{\alpha(m_h - M_0)} \left( \frac{t - t_h}{c} + 1\right)^{-p}}
#' Where \eqn{K, \alpha, c > 0}, and \eqn{p \geq 1} are the ETAS parameters, \eqn{t} (`t`) is the time
#' at which the function is evaluated, considering the past observation \eqn{t_h, m_h} (`th, mh`).
#' The function returns 0 if \eqn{t_h > t}.
#' If \eqn{t} is a scalar and \eqn{t_h, m_h} are vectors than the function returns a vector,
#' same if \eqn{t} is a vector and \eqn{t_h, m_h} are scalars, or if \eqn{t, t_h, m_h} are vectors of the same length.
#'
#' Do not use if \eqn{t} and \eqn{t_h, m_h} are vectors of different lengths.
#' @examples
gt <- function(theta, t, th, mh, M0){
  output <- rep(0,length(th))
  t.diff <- t - th
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    log.out <- log(theta$K) + theta$alpha*(mh[!neg] - M0)  - theta$p*log(1 + t.diff[!neg]/theta$c)
    output[!neg] <- exp(log.out)
  }
  else{
    output
  }
  output
}


#' ETAS conditional intensity - used by `inlabru`
#' @description Function to calculate the value of the ETAS model conditional intensity at a specified time given the history of the process.
#' @param theta ETAS parameters as `list` with names `mu`, `K`, `alpha`, `c`, `p`
#' @param t Time at which the conditional intensity is evaluated, `scalar`
#' @param th Time of the events in the history of the process, `vector`
#' @param mh Magnitudes of the events in the history of the process, `vector`
#' @param M0 Minimum magnitude threshold
#'
#' @return Value of the ETAS conditional intensity calculated at `t` with history `th, mh`, `scalar`
#' @export
#' @details
#' The function takes a single value `t` and returns the ETAS conditional intensity calculated at `t`
#' with history `th, mh`. The ETAS conditional intensity is given by
#' \deqn{\lambda(t | \mathcal H_t} = \mu + \sum_{h: (t_h,m_h) \in \mathcal H_t} K e^{\alpha(m_h - M_0)} \left( \frac{t - t_h}{c} + 1\right)^{-p}}
#'
#' Do not use when `t` is a vector.
#' @examples
cond.lambda <- function(theta, t, th, mh, M0){
  if(is.null(th) | all(th > t)){
    theta$mu
  }
  theta$mu + sum(gt(theta, t, th, mh, M0))
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
#' @param T2 End of temporal model domain.
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
  th + theta$c*( ( 1 - omega * (1/theta$c)*(theta$p - 1) )^( -1/(theta$p - 1) ) - 1)
}

#' Title
#'
#' @param th
#' @param t
#' @param ti
#' @param mi
#' @param M0
#'
#' @return
#' @export
#'
#' @examples
trigger <- function(th, t, ti, mi, M0){
  output <- rep(0,length(t))
  t.diff <- t - ti
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    log.out <- log(th[2]) + th[3]*(mi - M0) - th[5]*log(1 + t.diff[!neg]/th[4])
    output[!neg] <- exp(log.out)
  }
  else{
    output
  }
  output
}
