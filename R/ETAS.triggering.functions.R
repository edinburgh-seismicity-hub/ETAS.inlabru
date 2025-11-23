################################################################################
## Functions for ETAS conditional intensity and integrated conditional intensity
################################################################################


#' @title ETAS triggering function - used by `inlabru`
#' @description The function returns the value of the ETAS triggering function
#'   at a specified time `t` for the points in the history `th, mh`
#' @param theta ETAS parameters as `list` with names `K`, `alpha`, `c`, `p`
#' @param t Time at which the function is calculated, `scalar` or `vector`
#' @param th  Time of events in the history in `[days, months,...]`, `scalar` or
#'   `vector`
#' @param mh Magnitude of events in the history, `scalar` or `vector`
#' @param M0 Minimum magnitude threshold, `scalar`
#'
#' @return value of the ETAS triggering function evaluated at `t` with history
#'   `th`, `mh`.
#' @export
#' @details The ETAS triggering function to be evaluated is \deqn{g(t - t_h |
#'   m_h) = K e^{\alpha(m_h - M_0)} \left( \frac{t - t_h}{c} + 1\right)^{-p}}
#'   Where \eqn{K, \alpha, c > 0}, and \eqn{p \geq 1} are the ETAS parameters,
#'   \eqn{t} (`t`) is the time at which the function is evaluated, considering
#'   the past observation \eqn{t_h, m_h} (`th, mh`).
#' The function returns 0 if \eqn{t_h > t}. If \eqn{t} is a scalar and \eqn{t_h,
#' m_h} are vectors than the function returns a vector, same if \eqn{t} is a
#'   vector and \eqn{t_h, m_h} are scalars, or if \eqn{t, t_h, m_h} are vectors
#'   of the same length.
#'
#' Do not use if \eqn{t} and \eqn{t_h, m_h} are vectors of different lengths.
gt <- function(theta, t, th, mh, M0) {
  if (is.list(theta)) {
    mu <- theta$mu
    alpha <- theta$alpha
    K <- theta$K
    c <- theta$c
    p <- theta$p
  } else {
    mu <- theta[1]
    K <- theta[2]
    alpha <- theta[3]
    c <- theta[4]
    p <- theta[5]
  }

  output <- rep(0, length(th))
  t.diff <- t - th
  neg <- t.diff <= 0
  if (sum(!neg) > 0) {
    log.out <- log(K) + alpha * (mh - M0) - p * log(1 + t.diff[!neg] / c)
    # log.out <- log(theta[2]) + theta[3]*(mh[!neg] - M0)  - theta[5]*log(1 +
    # t.diff[!neg]/theta[4])
    output[!neg] <- exp(log.out)
  } else {
    output
  }
  output
}


#' @title ETAS conditional intensity - used by `inlabru`
#' @description Function to calculate the value of the ETAS model conditional
#'   intensity at a specified time given the history of the process.
#' @param theta ETAS parameters as `list` with names `mu`, `K`, `alpha`, `c`,
#'   `p`
#' @param t Time at which the conditional intensity is evaluated, `scalar`
#' @param th Time of the events in the history of the process, `vector`
#' @param mh Magnitudes of the events in the history of the process, `vector`
#' @param M0 Minimum magnitude threshold
#'
#' @return Value of the ETAS conditional intensity calculated at `t` with
#'   history `th, mh`, `scalar`
#' @export
#' @details
#' The function takes a single value `t` and returns the ETAS conditional
#' intensity calculated at `t` with history `th, mh`. The ETAS conditional
#' intensity is given by
#' \deqn{\lambda(t | \mathcal H_t) = \mu + \sum_{h: (t_h,m_h) \in \mathcal H_t}
#' K e^{\alpha(m_h - M_0)} \left( \frac{t - t_h}{c} + 1\right)^{-p}}
#'
#' Do not use when `t` is a vector.
cond_lambda <- function(theta, t, th, mh, M0) {
  if (is.null(th) || all(th > t)) {
    theta$mu
  }
  theta$mu + sum(gt(theta, t, th, mh, M0))
}


#' Logarithm of the integral of the ETAS triggering function
#'
#' @usage log_Lambda_h(theta, th, mh, M0, T1, T2)
#'
#' @param theta ETAS parameters `data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`.
#' @param th Time of parent event.
#' @param mh Magnitude of parent event
#' @param M0 Minimum magnitude threshold
#' @param T1 Start of temporal model domain.
#' @param T2 End of temporal model domain.
#'
#' @return Logarithm of the integral of the ETAS triggering function
#' @export log_Lambda_h
log_Lambda_h <- function(theta, th, mh, M0, T1, T2) {
  T.low <- pmax(T1, th)

  gamma.l <- (T.low - th) / theta$c
  gamma.u <- (T2 - th) / theta$c
  w.l <- (1 + gamma.l)^(1 - theta$p)
  w.u <- (1 + gamma.u)^(1 - theta$p)
  # output
  (log(theta$K) + theta$alpha * (mh - M0) + log(theta$c) - log(theta$p - 1) +
    log1p(w.l - 1) + log1p(-w.u / w.l))
}
