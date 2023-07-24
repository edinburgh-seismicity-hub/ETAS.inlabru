#' Retrieve posterior distribution of ETAS parameters
#'
#' @param input.list input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @return A `list` of two elements:
#' \itemize{
#' \item `post.df` : `data.frame` with the posterior distributions of the parameters with columns `x` (value of the parameter), `y` (value of the posterior), `param` (parameter name)
#' \item `post.plot` : `ggplot` object showing the posterior distribution of each parameter}
#' @export
get_posterior_param <- function(input.list) {
  post.mu <- data.frame(
    INLA::inla.tmarginal(
      input.list$link.functions$mu,
      input.list$model.fit$marginals.fixed$th.mu
    ),
    param = "mu"
  )
  post.K <- data.frame(
    INLA::inla.tmarginal(
      input.list$link.functions$K,
      input.list$model.fit$marginals.fixed$th.K
    ),
    param = "K"
  )
  post.alpha <- data.frame(
    INLA::inla.tmarginal(
      input.list$link.functions$alpha,
      input.list$model.fit$marginals.fixed$th.alpha
    ),
    param = "alpha"
  )
  post.c <- data.frame(
    INLA::inla.tmarginal(
      input.list$link.functions$c_,
      input.list$model.fit$marginals.fixed$th.c
    ),
    param = "c"
  )
  post.p <- data.frame(
    INLA::inla.tmarginal(
      input.list$link.functions$p,
      input.list$model.fit$marginals.fixed$th.p
    ),
    param = "p"
  )
  post.df <- rbind(post.mu, post.K, post.alpha, post.c, post.p)
  post.plot <- ggplot2::ggplot(post.df, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$param),
      scales = "free",
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::xlab("param") +
    ggplot2::ylab("pdf")
  list(
    post.df = post.df,
    post.plot = post.plot
  )
}


#' Sample from the posterior of the ETAS parameters
#'
#' @param input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp The number of samples to draw from the posteriors
#' @param max.batch Maximum number of posterior samples to be generated simultaneously. Default is 1000.
#' @param ncore Deprecated argument for controlling parallelism. Use
#' `future::plan(future::multisession, workers = ncore)` (or similar) to configure
#' parallelism in your code instead.
#'
#' @return `data.frame` of posterior samples with `nrow = n.samp` and columns `mu, K, alpha, c, p` corresponding to ETAS parameters.
#' @export
post_sampling <- function(input.list, n.samp, max.batch = 1000, ncore = NULL) {
  if (!is.null(ncore)) {
    lifecycle::deprecate_soft(
      "1.1.1.9001",
      "post_sampling(ncore)",
      I("future::plan(future::multisession, workers = ncore) in your code")
    )
  }
  if (n.samp > max.batch) {
    n.batch <- floor(n.samp / max.batch)
    if (n.samp - max.batch * n.batch == 0) {
      n.samp.per.batch <- rep(max.batch, n.batch)
    } else {
      n.samp.per.batch <- c(rep(max.batch, n.batch), n.samp - max.batch * n.batch)
    }
    post.samp.list <- future.apply::future_lapply(
      n.samp.per.batch,
      function(n.samp.batch) {
        inla.seed <- sample.int(n = .Machine$integer.max, size = 1)
        inlabru::generate(
          input.list$model.fit,
          data.frame(),
          ~ c(
            input.list$link.functions$mu(th.mu),
            input.list$link.functions$K(th.K),
            input.list$link.functions$alpha(th.alpha),
            input.list$link.functions$c_(th.c),
            input.list$link.functions$p(th.p)
          ),
          n.samples = n.samp.batch,
          seed = inla.seed
        )
      },
      future.seed = TRUE
    )
    post.samp <- do.call(cbind, post.samp.list)
  } else {
    post.samp <- inlabru::generate(
      input.list$model.fit,
      data.frame(),
      ~ c(
        input.list$link.functions$mu(th.mu),
        input.list$link.functions$K(th.K),
        input.list$link.functions$alpha(th.alpha),
        input.list$link.functions$c_(th.c),
        input.list$link.functions$p(th.p)
      ),
      n.samples = n.samp
    )
  }
  data.frame(
    mu = post.samp[1, ],
    K = post.samp[2, ],
    alpha = post.samp[3, ],
    c = post.samp[4, ],
    p = post.samp[5, ]
  )
}


#' Plot the posterior densities of the ETAS parameters
#'
#' @param input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp The number of samples to draw from the posteriors for the plot
#' @param post.samp `data.frame` with columns mu, K, alpha, c, p and rows corresponding to different posterior samples. When `NULL` the function samples the joint posterior distribution `n.samp` times. The default is `NULL`.
#' @param max.batch parameter of [post_sampling] function to be used in case `post.samp = NULL`
#' @return `list`: with elements
#' \itemize{
#' \item `post.samp.df`:`data.frame` of posterior samples with `nrow = n.samp` and columns `mu, K, alpha, c, p` corresponding to ETAS parameters. If `post.samp` is not `NULL` it returns `post.samp`
#' \item `pair.plot`: `ggplot` object reporting the pair plot between parameters samples. It is obtained using the `ggpairs` function of the `Ggally` library
#' }
#' @export
post_pairs_plot <- function(input.list = NULL, n.samp = NULL,
                            post.samp = NULL, max.batch = 1000) {
  if (is.null(input.list) & is.null(post.samp)) {
    stop("input.list and post.samp are missing, please provide at least one of the two")
  }
  if (is.null(post.samp)) {
    if (is.null(n.samp)) {
      stop("n.samp is missing, please provide a number of samples from the posterior")
    }
    post.samp <- post_sampling(input.list, n.samp = n.samp, max.batch = max.batch)
  }
  pair.plot <- GGally::ggpairs(post.samp,
    labeller = ggplot2::label_parsed,
    lower = list(continuous = "density")
  )
  return(list(
    post.samp.df = post.samp,
    pair.plot = pair.plot
  ))
}


#######
#' Calculate the integral of the ETAS conditional intensity
#'
#' @description Calculate the number of events in a time interval T1 to T2 given imposed events and ETAS
#'
#' @param th.mu Background rate, `mu`, on the internal parameter scale
#' @param th.K ETAS triggering parameter `K` on the internal parameter scale
#' @param th.alpha ETAS triggering parameter `alpha` on the internal parameter scale
#' @param th.c ETAS triggering parameter `c` on the internal parameter scale
#' @param th.p ETAS triggering parameter `p` on the internal parameter scale
#' @param T1 Start of temporal model domain.
#' @param T2 End of temporal model domain.
#' @param M0 Minimum magnitude threshold
#' @param Ht History of the process, or set of known events in the interval. It must be a `data.frame` with columns `ts` (time) and `magnitudes` (magnitudes).
#' @param link.functions `list` of functions to transform the parameters from the internal scale to the ETAS scale
#'
#' @return Integral of the ETAS conditional intensity between `T1` and `T2` with minimum magnitude `M0`,
#' i.e. the expected number of events.
#' @export
lambda_N <- function(th.mu, th.K, th.alpha, th.c, th.p, T1, T2, M0, Ht,
                     link.functions) {
  theta_etas <- list(
    mu = link.functions$mu(th.mu[1]),
    K = link.functions$K(th.K[1]),
    alpha = link.functions$alpha(th.alpha[1]),
    c = link.functions$c_(th.c[1]),
    p = link.functions$p(th.p[1])
  )


  return(theta_etas$mu * (T2 - T1) + sum(exp(log_Lambda_h(
    theta = theta_etas,
    th = Ht$ts,
    mh = Ht$magnitudes,
    M0 = M0,
    T1 = T1, T2 = T2
  ))))
}

#######
## MN: DESCRIPTION: Samples the distribution of distributions of N
## MN: Arguments: input.list combining input file with bru model output
## MN: Returns: (i) DataFrame with summary statistics for N
##              (ii) Plot of the expected number of events
##              (iii) Plot with shaded stack of all distributions of N and a normal about the mean model

#' Plot the posterior distribution of the expected number of events
#'
#' @param input.list Which has combined the input file (for link functions) and bru output (for marginals)
#' @param domain.extension Percentage of posterior quantiles to extend the domain specified as `scalar`. Default is set to 0.10.
#'
#' @return A `list` of three objects:
#' \itemize{
#' \item `post.df`: `data.frame` containing posterior informations on the posterior distribution of the number of events
#' \item `post.plot` : `ggplot` object showing the posterior distribution of the expected number of events
#' \item `post.plot.shaded` : `ggplot` object showing the posterior distribution of the expected number of events, shaded region corresponds to the 0.025 and 0.975 quantiles of the distribution of the distribution of the number of events.}
#' @export
get_posterior_N <- function(input.list, domain.extension = 0.10) {
  lambda.N.post <- predict(
    input.list$model.fit, # model fit
    data.frame(), # data (empty because the history of the process is passed to the function below directly)
    ~ lambda_N(
      th.mu, th.K, th.alpha, th.c, th.p,
      input.list$T12[1], input.list$T12[2], input.list$M0,
      input.list$catalog.bru,
      input.list$link.functions
    )
  ) # target function
  N.seq <- floor(lambda.N.post$q0.025 - lambda.N.post$q0.025 * domain.extension):ceiling(lambda.N.post$q0.975 + lambda.N.post$q0.975 * domain.extension)
  N.post.df <- predict(
    input.list$model.fit, data.frame(),
    ~ data.frame(
      N = N.seq,
      pdf = dpois(
        N.seq,
        lambda_N(
          th.mu = th.mu,
          th.K = th.K,
          th.alpha = th.alpha,
          th.c = th.c,
          th.p = th.p,
          T1 = input.list$T12[1],
          T2 = input.list$T12[2],
          M0 = input.list$M0,
          Ht = input.list$catalog.bru,
          link.functions = input.list$link.functions
        )
      )
    )
  )
  N.post.plot <- ggplot2::ggplot(N.post.df, ggplot2::aes(x = .data$N, y = .data$mean)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = nrow(input.list$catalog.bru), linetype = 2) +
    ggplot2::ylab("pdf")

  N.post.plot.shaded <- N.post.plot +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        xmax = .data$N,
        xmin = .data$N,
        ymin = .data$q0.025,
        ymax = .data$q0.975
      ),
      alpha = 0.2,
      fill = "blue"
    )
  list(
    post.df = N.post.df,
    post.plot = N.post.plot,
    post.plot.shaded = N.post.plot.shaded
  )
}
