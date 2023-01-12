





#' Function to plot the ETAS triggering function corresponding to different posterior samples
#'
#' @param list.input structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param magnitude Magnitude of the event for which the triggering function is calculated, `scalar` (`default = 4`).
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bund of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`)
#'
#' @return `ggplot` object with grey lines representing the triggering function for each posterior sample.
#' Black lines representing the 0.025 and 0.975 quantiles of the function values calculated for each posterior sample.
#' Horizontal red lines represents the 0.025 and 0.975 quantiles of the sampled background rates.
#' @export
#'
#' @examples
triggering_fun_plot <- function(list.input, magnitude = 4, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  post.samp <- generate(list.input$model.fit, data.frame(),
                        ~ c(list.input$link.functions$mu(th.mu),
                            list.input$link.functions$K(th.K),
                            list.input$link.functions$alpha(th.alpha),
                            list.input$link.functions$c_(th.c),
                            list.input$link.functions$p(th.p)
                        ), n.samples = n.samp)
  post.samp <- t(post.samp)

  trig.eval <- lapply(1:nrow(post.samp),
                      \(x) gt(theta = post.samp[x,],
                              t = t.eval,
                              th = 0,
                              mh = magnitude,
                              M0 = list.input$M0))
  trig.cols <- as.matrix(bind_cols(trig.eval))
  trig.lower.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.025)))
  trig.upper.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.975)))
  mu.lower.quant <- quantile(post.samp[,1], 0.025)
  mu.upper.quant <- quantile(post.samp[,1], 0.975)
  #omori.eval <- bind_rows(omori.eval)
  output.plot <- ggplot()
  for(i in 1:ncol(trig.cols)){
    trig.eval.i <- trig.cols[,i]
    df.trig <- data.frame(t = t.eval,
                          trig = trig.eval.i)
    output.plot <- output.plot +
      geom_line(data = df.trig,
                aes(x = t, y = trig),
                color= 'grey', alpha = 0.5)
  }
  output.plot +
    geom_line(aes(x = t.eval, y = trig.lower.quant)) +
    geom_line(aes(x = t.eval, y = trig.upper.quant)) +
    geom_hline(yintercept = mu.lower.quant, color = 'red') +
    geom_hline(yintercept = mu.upper.quant, color = 'red') +
    theme_bw() +
    xlab("Time") +
    ylab("Event rate")
}


#' Function to plot the ETAS triggering function corresponding to different posterior samples
#'
#' @param list.input structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param magnitude Magnitude of the event for which the triggering function is calculated, `scalar` (`default = 4`).
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bund of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`)
#'
#' @return `ggplot` object with grey lines representing the triggering function for each posterior sample.
#' Black lines representing the 0.025 and 0.975 quantiles of the function values calculated for each posterior sample.
#' Horizontal red lines represents the 0.025 and 0.975 quantiles of the sampled background rates.
#' @export
#'
#' @examples
triggering_fun_plot_priors <- function(list.input, magnitude = 4, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  prior.samp <- cbind( list.input$link.functions$mu(rnorm(n.samp)),
                      list.input$link.functions$K(rnorm(n.samp)),
                      list.input$link.functions$alpha(rnorm(n.samp)),
                      list.input$link.functions$c(rnorm(n.samp)),
                      list.input$link.functions$p(rnorm(n.samp)))

  trig.eval <- lapply(1:nrow(prior.samp),
                      \(x) gt(theta = prior.samp[x,],
                              t = t.eval,
                              th = 0,
                              mh = magnitude,
                              M0 = list.input$M0))
  trig.cols <- as.matrix(bind_cols(trig.eval))
  trig.lower.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.025)))
  trig.upper.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.975)))
  mu.lower.quant <- quantile(prior.samp[,1], 0.025)
  mu.upper.quant <- quantile(prior.samp[,1], 0.975)
  #omori.eval <- bind_rows(omori.eval)
  output.plot <- ggplot()
  for(i in 1:ncol(trig.cols)){
    trig.eval.i <- trig.cols[,i]
    df.trig <- data.frame(t = t.eval,
                          trig = trig.eval.i)
    output.plot <- output.plot +
      geom_line(data = df.trig,
                aes(x = t, y = trig),
                color= 'grey', alpha = 0.5)
  }
  output.plot +
    geom_line(aes(x = t.eval, y = trig.lower.quant)) +
    geom_line(aes(x = t.eval, y = trig.upper.quant)) +
    geom_hline(yintercept = mu.lower.quant, color = 'red') +
    geom_hline(yintercept = mu.upper.quant, color = 'red') +
    theme_bw() +
    xlab("Time") +
    ylab("Event rate")
}



#' Function to calculate the omori's law
#'
#' @param theta ETAS parameters (`list(mu = mu, K = K, alpha = alpha, c = c, p = p`), only parameters `c` and `p` are used
#' @param t Time at which the Omori's law is evaluated
#' @param ti Time of the event in the history
#'
#' @return Value of the Omori's law at point `t` for and event happened in `ti`
#'
#' @examples
omori <- function(theta, t, ti){
  output <- rep(0,length(t))
  t.diff <- t - ti
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    log.out <- - theta[5]*log(1 + t.diff[!neg]/theta[4])
    output[!neg] <- exp(log.out)
  }
  else{
    output
  }
  output
}

#' Function to plot the Omori's law corresponding to different posterior samples
#'
#' @param list.input structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bund of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`).
#'
#' @return
#' @export
#'
#' @examples
omori_plot <- function(list.input, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  prior.samp <- cbind( list.input$link.functions$mu(rnorm(n.samp)),
                       list.input$link.functions$K(rnorm(n.samp)),
                       list.input$link.functions$alpha(rnorm(n.samp)),
                       list.input$link.functions$c(rnorm(n.samp)),
                       list.input$link.functions$p(rnorm(n.samp)))

  omori.eval <- lapply(1:nrow(post.samp),
                       \(x) omori(th = post.samp[x,],
                                  t = t.eval,
                                  ti = 0))
  omori.cols <- as.matrix(bind_cols(omori.eval))
  omori.lower.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.025)))
  omori.upper.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.975)))

  #omori.eval <- bind_rows(omori.eval)
  output.plot <- ggplot()
  for(i in 1:ncol(omori.cols)){
    omori.eval.i <- omori.cols[,i]
    df.omori <- data.frame(t = t.eval,
                           omori = omori.eval.i)
    output.plot <- output.plot +
      geom_line(data = df.omori,
                aes(x = t, y = omori),
                color= 'grey', alpha = 0.5)
  }
  output.plot +
    geom_line(aes(x = t.eval, y = omori.lower.quant)) +
    geom_line(aes(x = t.eval, y = omori.upper.quant)) +
    theme_bw() +
    xlab("Time") +
    ylab("Event rate")
}

#' Function to plot the Omori's law corresponding to different posterior samples
#'
#' @param list.input structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bund of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`).
#'
#' @return
#' @export
#'
#' @examples
omori_plot_priors <- function(list.input, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  post.samp <- generate(list.input$model.fit, data.frame(),
                        ~ c(list.input$link.functions$mu(th.mu),
                            list.input$link.functions$K(th.K),
                            list.input$link.functions$alpha(th.alpha),
                            list.input$link.functions$c_(th.c),
                            list.input$link.functions$p(th.p)
                        ), n.samples = n.samp)
  post.samp <- t(post.samp)

  omori.eval <- lapply(1:nrow(post.samp),
                       \(x) omori(th = post.samp[x,],
                                  t = t.eval,
                                  ti = 0))
  omori.cols <- as.matrix(bind_cols(omori.eval))
  omori.lower.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.025)))
  omori.upper.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.975)))

  #omori.eval <- bind_rows(omori.eval)
  output.plot <- ggplot()
  for(i in 1:ncol(omori.cols)){
    omori.eval.i <- omori.cols[,i]
    df.omori <- data.frame(t = t.eval,
                           omori = omori.eval.i)
    output.plot <- output.plot +
      geom_line(data = df.omori,
                aes(x = t, y = omori),
                color= 'grey', alpha = 0.5)
  }
  output.plot +
    geom_line(aes(x = t.eval, y = omori.lower.quant)) +
    geom_line(aes(x = t.eval, y = omori.upper.quant)) +
    theme_bw() +
    xlab("Time") +
    ylab("Event rate")
}
