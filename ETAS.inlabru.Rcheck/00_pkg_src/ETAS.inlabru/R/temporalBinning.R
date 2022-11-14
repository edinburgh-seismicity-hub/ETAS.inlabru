# find breaks point for grid
#' Title
#'
#' @param tt_ List of the event times a grid is needed for [days].
#' @param T2_ End of temporal domain [days].
#' @param coef_ TimeBinning parameter:
#' @param delta_ TimeBinning parameter:
#' @param N_exp_ TimeBinning parameter: Number of bins
#'
#' @return
#' @export
#'
#' @examples
breaks_exp <- function(tt_, T2_, coef_ = 2, delta_, N_exp_ = 10){

  # Generate the breaks for each event
  tt_breaks <- tt_ + delta_*((1 + coef_)^(0:N_exp_))

  # Remove the breaks ending before the end of the temporal domain, T2.
  tt_breaks <- tt_breaks[tt_breaks < T2_]

  # If domain smaller than min threshold, only return bounds
  if(T2_ - tt_ < delta_){
    return(c(tt_, T2_))
  }
  if(T2_ - tt_breaks[length(tt_breaks)] < delta_){
    tt_breaks[length(tt_breaks)] = T2_
  }
  if(tt_breaks[length(tt_breaks)] < T2_){
    tt_breaks <- c(tt_breaks, T2_)
  }
  return(c(tt_,tt_breaks))
}


#' Generate a set of time bins for a specific event and return.
#'
#' @param data.point
#' @param coef.t TimeBinning parameter:
#' @param delta.t TimeBinning parameter:
#' @param T2. End of the temporal domain [days].
#' @param displaygrid Boolean variable - whether to plot the grid
#' @param N.exp. TimeBinning parameter:
#'
#' @return A set of time bins aggregated over all events MORE DETAIL
#' @export
#'
#' @examples
#' ## EXAMPLE 1
#' events <- data.frame( ts=c(0,1 , 3 ), idx.p=c(1,2,3) )
#' T2 <- 20
#' N.exp <- 8
#' delta.t <- 1
#' coef.t <- 1
#' time.grid(events, coef.t, delta.t, T2, displaygrid = FALSE, N.exp)
time.grid <- function(data.point, coef.t, delta.t,
                      T2., displaygrid = FALSE, N.exp.){

  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # spatial bins
  if(displaygrid){
    Plot_grid(xx = xx., yy = yy., delta_ = delta., n.layer = n.layer.,
              bdy_ =  bdy., min.edge = min.edge.)
  }

  # time bins
  # find bins break points
  t_b <- breaks_exp(tt., T2., coef_ = coef.t, delta_ = delta.t, N_exp_ = N.exp.)
  time.bins <- data.frame(t.start = t_b[-length(t_b)],
                          t.end = t_b[-1]) %>%
    mutate(t.bin.name = paste0(round(t.start,3),'-',round(t.end,3)))
  if(nrow(time.bins) - 1 == 0){
    time.bins$t.ref_layer = paste0('last-',idx.p)
  }
  else{
    time.bins$t.ref_layer = c(1:(nrow(time.bins) - 1), paste0('last-',idx.p))
  }
  return( cbind(time.bins, data.point, row.names = NULL) )
}

#' Title
#'
#' @param param_
#' @param time.df
#'
#' @return
#' @export
#'
#' @examples
It_df <- function(param_, time.df){
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_c <- param_[4]
  param_p <- param_[5]
  T.l <- pmax(tth, T1b) #sapply(1:length(tth), \(x) max(tth[x], T1b[x]))
  fun.l <- (1 + (T.l - tth)/param_c)^(1-param_p)
  fun.u <- (1 + (T2b - tth)/param_c)^(1-param_p)
  ( param_c/ (param_p - 1) )* ( fun.l - fun.u )
}


#' Title
#'
#' @param param.
#' @param list.input_
#'
#' @return
#' @export
#'
#' @examples
compute.grid <- function(param., list.input_){

  It.vec <- It_df(param_ = param., time.df = list.input_$time.sel)

  It.vec[list.input_$Imapping]
}

