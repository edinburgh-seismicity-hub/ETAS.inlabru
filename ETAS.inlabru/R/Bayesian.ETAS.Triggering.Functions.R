#' Conditional intensity - used by bayesianETAS
#' Used for comparing output of inlabru with Bayesian ETAS
#'
#' @param th
#' @param t
#' @param ti.v
#' @param mi.v
#' @param M0 Minimum magnitude threshold
#'
#' @return
#' @export
#'
#' @examples
lambda_be <- function(th, t, ti.v, mi.v, M0){
  if(is.null(ti.v) | all(ti.v > t)){
    th[1]
  }
  th[1] + sum(gt(th, t, ti.v, mi.v, M0))
}

#' Time triggering function - used by bayesianETAS
#' Used for comparing output of inlabru with Bayesian ETAS
#' MN: TODO Cross-reference the paper
#'
#' @param th
#' @param t
#' @param ti
#' @param mi
#' @param M0 Minimum magnitude threshold
#'
#' @return
#' @export
#'
#' @examples
gt_be <- function(th, t, ti, mi, M0){
  # set 0
  output <- rep(0,length(ti))
  before <- ti < t
  # update only ti < t
  if(sum(before) > 0){
    log.out <- log(th[2]) + th[3]*(mi[before] - M0)  + log(th[5] - 1) + (th[5] - 1)*log(th[4]) - th[5]*log(t - ti[before] + th[4])
    output[before] <- exp(log.out)
  }
  else{
    output
  }
  output
}

#' integrated triggering function - used by bayesianETAS
#' Used for comparing output of inlabru with Bayesian ETAS
#'
#' @param th
#' @param ti
#' @param mi
#' @param M0 Minimum magnitude threshold
#' @param T1 Start of temporal model domain.
#' @param T2 End of temporal model domain.
#'
#' @return
#' @export
#'
#' @examples
log.Lambda.h_be <- function(th, ti, mi, M0, T1, T2){
  th <- as.numeric(th)
  T.low <- pmax(ti, T1)#sapply(ti, \(x) max(T1, x))
  log(th[2]) + th[3]*(mi - M0) + log((th[4]/(T.low - ti + th[4]))^(th[5] - 1) - (th[4]/(T2 - ti + th[4]))^(th[5] - 1))
}
