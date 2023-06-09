
#############################
#### Injection rate function
#' @title Injection Rate function
#'
#' @description
#' Forward time integrated function for exponential rate decay
#'
#' @param a Event rate per unit volume injected
#' @param V.i Injected volume
#' @param tau Decau rate `[days]`
#' @param T.i Time of injection event
#' @param T2 TODO
#'
#' @return
#'
#' @examples
IntInjecIntensity <- function(a=50, V.i=1, tau=10, T.i, T2){
  expected.injection.events <- - tau*V.i*a* ( exp(-(T2-T.i)/tau ) -1 )
  return(expected.injection.events)
}

#' @title Returns end time given a ...
#'
#' @description
#' Returns end time given a ...
#'
#' @param a Event rate per unit volume injected
#' @param V.i Injected volume
#' @param tau Decay rate `[days]`
#' @param T.i Time of injection event
#' @param number.injected.events TODO
#'
#' @return
#' @export
#'
#' @examples
Inv.IntInjecIntensity <- function(a=50, V.i=1, tau=10, T.i, number.injected.events){
  endTime <- T.i - tau*log(1 - number.injected.events / (tau*V.i*a))
  return(endTime)
}

#' Title
#'
#' @param a Induced event rate per unit volume.
#' @param V.i Injected volume
#' @param tau Decay rate `[days]`.
#' @param beta.p Related to the b-value via `b ln(10)`.
#' @param M0 Minimum magnitude threshold.
#' @param T.i Time of injection `[days]`.
#' @param T2 End of temporal model domain `[days]`.
#'
#' @return Catalogue of parent events induced by injection data.frame(times, magnitudes)
#' @export
#'
#' @examples
sample.temoral.injection.events <- function(a=50, V.i=1, tau=10, beta.p, M0, T.i, T2){
  bound.l <- 0 #It(th.p, th, T)
  bound.u <- IntInjecIntensity(a=a, V.i=V.i, tau=tau, T.i=T.i, T2=T2)
  n.ev <- rpois( 1, bound.u  )
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  sample.ts <- Inv.IntInjecIntensity(a=a, V.i=V.i, tau=tau, T.i=T.i, number.injected.events=unif.s)

  samp.mags <- rexp(n.ev, rate = beta.p) + M0

  samp.points <- data.frame(ts = sample.ts, magnitudes = samp.mags)
  return(samp.points[!is.na(samp.points$ts),])
}

