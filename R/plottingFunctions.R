



#' Function to plot the ETAS triggering function corresponding to different posterior samples
#'
#' @param input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param post.samp `data.frame` containing posterior samples of the parameters. If `NULL`, then `n.samp` samples are
#' generated. If `n.samp` is different from `nrow(post.samp)` then `n.samp` rows are uniformly sampled from `post.samp`. Default is `NULL`
#' @param magnitude Magnitude of the event for which the triggering function is calculated, `scalar` (`default = 4`).
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bound of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`)
#'
#' @return `ggplot` object with grey lines representing the triggering function for each posterior sample.
#' Black lines representing the 0.025 and 0.975 quantiles of the function values calculated for each posterior sample.
#' Horizontal red lines represents the 0.025 and 0.975 quantiles of the sampled background rates.
#' @export
triggering_fun_plot <- function(input.list, post.samp = NULL, n.samp = 10,
                                magnitude = 4, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  if(is.null(post.samp)){
    post.samp <- inlabru::generate(
      input.list$model.fit,
      data.frame(),
      ~ c(input.list$link.functions$mu(th.mu),
          input.list$link.functions$K(th.K),
          input.list$link.functions$alpha(th.alpha),
          input.list$link.functions$c_(th.c),
          input.list$link.functions$p(th.p)
      ),
      n.samples = n.samp)
    post.samp <- t(post.samp)
  }
  if((!is.null(post.samp))&(!is.null(n.samp))){
    if(nrow(post.samp) != n.samp){
      idx.samp <- sample(1:nrow(post.samp), n.samp, replace = TRUE)
      post.samp <- post.samp[idx.samp,]
    }
  }
  trig.eval <- lapply(seq_len(nrow(post.samp)),
                      \(x) gt(theta = post.samp[x,],
                              t = t.eval,
                              th = 0,
                              mh = magnitude,
                              M0 = input.list$M0))
  trig.cols <- do.call(cbind, trig.eval)
  trig.lower.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.025)))
  trig.upper.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.975)))
  mu.lower.quant <- quantile(post.samp[,1], 0.025)
  mu.upper.quant <- quantile(post.samp[,1], 0.975)

  trig.eval <- lapply(
    seq_along(trig.eval),
    function(k) {
      data.frame(
        sample = k,
        t = t.eval,
        trig = trig.eval[[k]]
      )
    }
  )
  trig.eval <- dplyr::bind_rows(trig.eval)

  output.plot <- ggplot2::ggplot()
  output.plot <- output.plot +
    ggplot2::geom_line(
      data = trig.eval,
      ggplot2::aes(
        x = .data$t,
        y = .data$trig,
        group = factor(.data$sample)
      ),
      color = 'grey',
      alpha = 0.5
    )

  output.plot +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = trig.lower.quant)) +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = trig.upper.quant)) +
    ggplot2::geom_hline(yintercept = mu.lower.quant, color = 'red') +
    ggplot2::geom_hline(yintercept = mu.upper.quant, color = 'red') +
    ggplot2::theme_bw() +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Event rate")
}


#' Function to plot the ETAS triggering function corresponding to different prior samples
#'
#' @param input.list structured input `list` with at least one element:
#' \itemize{
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param magnitude Magnitude of the event for which the triggering function is calculated, `scalar` (`default = 4`).
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bound of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`)
#'
#' @return `ggplot` object with grey lines representing the triggering function for each posterior sample.
#' Black lines representing the 0.025 and 0.975 quantiles of the function values calculated for each posterior sample.
#' Horizontal red lines represents the 0.025 and 0.975 quantiles of the sampled background rates.
#' @export
triggering_fun_plot_prior <- function(input.list, magnitude = 4, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  prior.samp <- cbind( input.list$link.functions$mu(rnorm(n.samp)),
                      input.list$link.functions$K(rnorm(n.samp)),
                      input.list$link.functions$alpha(rnorm(n.samp)),
                      input.list$link.functions$c_(rnorm(n.samp)),
                      input.list$link.functions$p(rnorm(n.samp)))

  trig.eval <- lapply(seq_len(nrow(prior.samp)),
                      \(x) gt(theta = prior.samp[x,],
                              t = t.eval,
                              th = 0,
                              mh = magnitude,
                              M0 = input.list$M0))
  trig.cols <- do.call(cbind, trig.eval)
  trig.lower.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.025)))
  trig.upper.quant <- apply(trig.cols, 1, \(x) quantile(x, c(0.975)))
  mu.lower.quant <- quantile(prior.samp[,1], 0.025)
  mu.upper.quant <- quantile(prior.samp[,1], 0.975)

  trig.eval <- lapply(
    seq_along(trig.eval),
    function(k) {
      data.frame(
        sample = k,
        t = t.eval,
        trig = trig.eval[[k]]
      )
    }
  )
  trig.eval <- dplyr::bind_rows(trig.eval)

  output.plot <- ggplot2::ggplot()
  output.plot <- output.plot +
    ggplot2::geom_line(
      data = trig.eval,
      ggplot2::aes(
        x = .data$t,
        y = .data$trig,
        group = factor(.data$sample)
      ),
      color = 'grey',
      alpha = 0.5
    )

  output.plot +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = trig.lower.quant)) +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = trig.upper.quant)) +
    ggplot2::geom_hline(yintercept = mu.lower.quant, color = 'red') +
    ggplot2::geom_hline(yintercept = mu.upper.quant, color = 'red') +
    ggplot2::theme_bw() +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Event rate")
}



#' Function to calculate Omori's law
#'
#' @param theta ETAS parameters (`list(mu = mu, K = K, alpha = alpha, c = c, p = p`), only parameters `c` and `p` are used
#' @param t Time at which Omori's law is evaluated
#' @param ti Time of the event in the history
#'
#' @return Value of Omori's law at point `t` for and event happened in `ti`
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

#' Function to plot Omori's law corresponding to different prior samples
#'
#' @param input.list structured input `list` with at least one element:
#' \itemize{
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bound of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`).
#'
#' @return A ggplot object
#' @seealso [create.input.list.temporal.noCatalogue()], [create.input.list.temporal.withCatalogue()]
#' @export
omori_plot_prior <- function(input.list, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  prior.samp <- cbind( input.list$link.functions$mu(rnorm(n.samp)),
                       input.list$link.functions$K(rnorm(n.samp)),
                       input.list$link.functions$alpha(rnorm(n.samp)),
                       input.list$link.functions$c_(rnorm(n.samp)),
                       input.list$link.functions$p(rnorm(n.samp)))

  omori.eval <- lapply(seq_len(nrow(prior.samp)),
                       \(x) omori(theta = prior.samp[x,],
                                  t = t.eval,
                                  ti = 0))
  omori.cols <- do.call(cbind, omori.eval)
  omori.lower.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.025)))
  omori.upper.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.975)))

  omori.eval <- lapply(
    seq_along(omori.eval),
    function(k) {
      data.frame(
        sample = k,
        t = t.eval,
        omori = omori.eval[[k]]
      )
    }
  )
  omori.eval <- dplyr::bind_rows(omori.eval)

  output.plot <- ggplot2::ggplot()
  output.plot <- output.plot +
    ggplot2::geom_line(
      data = omori.eval,
      ggplot2::aes(
        x = .data$t,
        y = .data$omori,
        group = factor(.data$sample)
      ),
      color = 'grey',
      alpha = 0.5
    )

  output.plot +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = omori.lower.quant)) +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = omori.upper.quant)) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Event rate")
}

#' Function to plot Omori's law corresponding to different posterior samples
#'
#' @param input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp Number of posterior samples, `integer` (`default = 10`).
#' @param t.end Upper bound of the x-axis, `scalar` (`default = 1`).
#' @param n.breaks Number of points between 0 and `t.end` to calculate the function, `integer` (`default = 100`).
#'
#' @return A ggplot object
#' @seealso [create.input.list.temporal.noCatalogue()], [create.input.list.temporal.withCatalogue()]
#' @export
omori_plot_posterior <- function(input.list, n.samp = 10, t.end = 1, n.breaks = 100){
  t.eval <- seq(1e-6, t.end, length.out = n.breaks)
  post.samp <- inlabru::generate(
    input.list$model.fit,
    data.frame(),
    ~ c(input.list$link.functions$mu(th.mu),
        input.list$link.functions$K(th.K),
        input.list$link.functions$alpha(th.alpha),
        input.list$link.functions$c_(th.c),
        input.list$link.functions$p(th.p)
    ),
    n.samples = n.samp)
  post.samp <- t(post.samp)

  omori.eval <- lapply(seq_len(nrow(post.samp)),
                       \(x) omori(theta = post.samp[x,],
                                  t = t.eval,
                                  ti = 0))
  omori.cols <- do.call(cbind, omori.eval)
  omori.lower.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.025)))
  omori.upper.quant <- apply(omori.cols, 1, \(x) quantile(x, c(0.975)))

  omori.eval <- lapply(
    seq_along(omori.eval),
    function(k) {
      data.frame(
        sample = k,
        t = t.eval,
        omori = omori.eval[[k]]
      )
    }
  )
  omori.eval <- dplyr::bind_rows(omori.eval)

  output.plot <- ggplot2::ggplot()
  output.plot <- output.plot +
    ggplot2::geom_line(data = omori.eval,
                       ggplot2::aes(x = t, y = omori, group = factor(sample)),
                       color= 'grey', alpha = 0.5)
  output.plot +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = omori.lower.quant)) +
    ggplot2::geom_line(ggplot2::aes(x = t.eval, y = omori.upper.quant)) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Event rate")
}
