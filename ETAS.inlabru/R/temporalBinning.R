
#' Find breaks point for 1D grid
#' @description `breaks_exp` return the breaks points of a one dimensional grid depending on three parameters, see details
#' @param start.grid Starting point of the grid, `scalar`.
#' @param end.grid End point of the grid, `scalar`.
#' @param coef.t TimeBinning parameter: \eqn{\delta > 0} determines the relative length of subsequent intervals, `scalar`.
#' @param delta.t TimeBinning parameter: \eqn{\Delta > 0} determines the length of the intervals, `scalar`.
#' @param N.exp. TimeBinning parameter: \eqn{n_{max} >0} determines the maximum number of intervals, `scalar`
#'
#' @details The grid is calculated as follows
#' \deqn{t, t + \Delta, t + \Delta(1 + \delta), t + \Delta(1 + \delta)^2,...., T}
#' where \eqn{t} is the `start.grid` argument, \eqn{T} is the `end.grid` argument, and \eqn{n_{max}} is the maximum value of the exponent
#' @return `vector` containing the grid points
#'
#' @examples
#' breaks_exp(start.grid = 1, end.grid = 100, coef.t = 1, delta.t = 1, N.exp. = 3)
#' breaks_exp(start.grid = 1, end.grid = 100, coef.t = 1, delta.t = 1, N.exp. = 10)
#'
#' @export
breaks_exp <- function(start.grid, end.grid, coef.t = 2, delta.t, N.exp. = 10){

  # Generate the breaks for each event
  tt_breaks <- start.grid + delta.t*((1 + coef.t)^(0:N.exp.))

  # Remove the breaks ending before the end of the temporal domain, T2.
  tt_breaks <- tt_breaks[tt_breaks < end.grid]

  # If domain smaller than min threshold, only return bounds
  if(end.grid - start.grid < delta.t){
    return(c(start.grid, end.grid))
  }
  if(end.grid - tt_breaks[length(tt_breaks)] < delta.t){
    tt_breaks[length(tt_breaks)] = end.grid
  }
  if(tt_breaks[length(tt_breaks)] < end.grid){
    tt_breaks <- c(tt_breaks, end.grid)
  }
  return(c(start.grid,tt_breaks))
}


#' Generate a set of time bins for a specific event.
#' @param data.point Point for which the binning is calculated, `list` with elements time (`ts, scalar`), event index (`idx.p, scalar`). Names are mandatory and should not be changed
#' @param coef.t TimeBinning parameter: look [breaks_exp()]
#' @param delta.t TimeBinning parameter: look [breaks_exp()]
#' @param N.exp. TimeBinning parameter: look [breaks_exp()]
#' @param T1. Start of the temporal domain, `scalar`
#' @param T2. End of the temporal domain `scalar`.
#'
#' @return
#' A `data.frame` with as many rows as the number of bins and fixed number of columns.
#' The columns are
#' \itemize{
#' \item `t.start` : starting point of the bin (minimum = `T1.`).
#' \item `t.end` : end point of the bin. (maximum = `T2.`).
#' \item `t.bin.name` : unique bin identifier.
#' \item `t.ref_layer` : bin identifier for calculations
#' \item `ts` : time provided in `data.point`
#' \item `idx.p` : identifier provided in `data.point`
#' }
#' The bins are only between `T1.` and `T2.` or containing `T1.`
#' @export time.grid
#'
#' @examples
#' ## EXAMPLE 1
#' event <- list( ts= 0, idx.p= 1 )
#' time.grid(data.point = event, coef.t = 1, delta.t = 0.1, N.exp. = 8, T1 = 1, T2 = 20)
time.grid <- function(data.point, coef.t, delta.t, N.exp.,
                      T1., T2.){
  # extract point information
  tt. <- data.point$ts
  idx.p <- data.point$idx.p

  # time bins
  # find bins break points
  t_b <- breaks_exp(start.grid = tt.,
                    end.grid = T2.,
                    coef.t = coef.t,
                    delta.t = delta.t,
                    N.exp. = N.exp.)
  # create data.frame with t.start (bin starting time), t.end (bin end time) and t.bin.name (unique name for each bin)
  time.bins <- data.frame(t.start = t_b[-length(t_b)],
                          t.end = t_b[-1]) %>%
    dplyr::mutate(t.bin.name = paste0(round(t.start,3),'-',round(t.end,3)))
  # if there is only 1 bin than it is the last
  if(nrow(time.bins) - 1 == 0){
    time.bins$t.ref_layer = paste0('last-',idx.p)
  }
  # assign a number to each bin except the last one.
  else{
    time.bins$t.ref_layer = c(1:(nrow(time.bins) - 1), paste0('last-',idx.p))
  }
  # remove bins before T1
  bin.before <- time.bins$t.end < T1.
  time.bins <- time.bins[!bin.before,]
  # change name to bins containing T1
  if(sum(time.bins$t.start <= T1. & time.bins$t.end >= T1.) > 0){
    time.bins$t.ref_layer[time.bins$t.start <= T1. & time.bins$t.end >= T1.] <- paste0('between-',idx.p)
    time.bins$t.start[time.bins$t.start <= T1. & time.bins$t.end >= T1.] <- T1.
  }
  # merge with data.frame taken in input
  cbind(time.bins, data.point, row.names = NULL)
}

#' Function to calculate the integral of the Omori's law
#'
#' @param param_ ETAS parameters vector (\eqn{\mu, K, \alpha, c, p}), only \eqn{c, p} are used.
#' @param time.df output of the function [time.grid()]
#'
#' @return `vector` of integral values between each bin provided in `time.df`
It_df <- function(param_, time.df){
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_c <- param_[4]
  param_p <- param_[5]
  T.l <- pmax(tth, T1b)
  fun.l <- (1 + (T.l - tth)/param_c)^(1-param_p)
  fun.u <- (1 + (T2b - tth)/param_c)^(1-param_p)
  ( param_c/ (param_p - 1) )* ( fun.l - fun.u )
}


#' Function to compute the integral of Omori's law efficiently
#'
#' @param param. ETAS parameters vector (\eqn{\mu, K, \alpha, c, p}), only \eqn{c, p} are used.
#' @param list.input_ `list` containing information to calculate the integrals efficiently. The list is created inside the `Temporal.ETAS` function
#' Two elements are used
#' \itemize{
#' \item `time.sel` selection of rows of the output of `time.grid` with unique `t.ref_layer` value, `data.frame`.
#' \item `Imapping` mapper between the unique names provided in `time.sel` and original rows of the output of `time.grid()`, `vector`.
#' }
#' @return `vector` with same length as `list.input_$Imapping` with the integral of the Omori's law for each bin.
compute.grid <- function(param., list.input_){

  It.vec <- It_df(param_ = param., time.df = list.input_$time.sel)

  It.vec[list.input_$Imapping]
}

