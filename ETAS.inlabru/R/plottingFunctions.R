#' Title
#'
#' @param xx
#' @param yy
#' @param delta_
#' @param n.layer
#' @param bdy_
#' @param min.edge
#'
#' @return
#' @export
#'
#' @examples
Plot_grid <- function (xx = xx., yy = yy., delta_ = delta., n.layer = n.layer.,
          bdy_ =  bdy., min.edge = min.edge.) {
  return( NULL )
}

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
                      \(x) trigger(th = post.samp[x,],
                                   t = t.eval,
                                   ti = 0,
                                   mi = magnitude,
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


#' Title
#'
#' @param th
#' @param t
#' @param ti
#'
#' @return
#' @export
#'
#' @examples
omori <- function(th, t, ti){
  output <- rep(0,length(t))
  t.diff <- t - ti
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    log.out <- - th[5]*log(1 + t.diff[!neg]/th[4])
    output[!neg] <- exp(log.out)
  }
  else{
    output
  }
  output
}

#' Plot the Omori function using samples from the posteriors containted in `list.input``
#'
#' @param list.input
#' @param n.samp
#' @param t.end
#' @param n.breaks
#'
#' @return
#' @export
#'
#' @examples
omori_plot <- function(list.input, n.samp = 10, t.end = 1, n.breaks = 100){
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
