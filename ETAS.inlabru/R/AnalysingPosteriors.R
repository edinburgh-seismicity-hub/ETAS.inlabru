library(GGally)
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
#'
#' @examples
get_posterior_param <- function(input.list){
  post.mu <- data.frame(inla.tmarginal(input.list$link.functions$mu,
                                       input.list$model.fit$marginals.fixed$th.mu),
                        param = 'mu')
  post.K <- data.frame(inla.tmarginal(input.list$link.functions$K,
                                      input.list$model.fit$marginals.fixed$th.K),
                       param = 'K')
  post.alpha <- data.frame(inla.tmarginal(input.list$link.functions$alpha,
                                          input.list$model.fit$marginals.fixed$th.alpha),
                           param = 'alpha')
  post.c <- data.frame(inla.tmarginal(input.list$link.functions$c_,
                                      input.list$model.fit$marginals.fixed$th.c),
                       param = 'c')
  post.p <- data.frame(inla.tmarginal(input.list$link.functions$p,
                                      input.list$model.fit$marginals.fixed$th.p),
                       param = 'p')
  post.df <- rbind(post.mu, post.K, post.alpha, post.c, post.p)
  post.plot <- ggplot(post.df, aes(x,y)) +
    geom_line() +
    facet_wrap(facets = vars(param), scales = 'free', labeller = label_parsed)+
    xlab('param') +
    ylab('pdf')
  list(post.df = post.df,
       post.plot = post.plot)
}


#' Sample from the posterior of the ETAS parameters
#'
#' @param input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp The number of samples to draw from the posteriors
#'
#' @return `data.frame` of posterior samples with `nrow = n.samp` and columns `mu, K, alpha, c, p` corresponding to ETAS parameters.
#' @export
#'
#' @examples
post_sampling <- function(input.list, n.samp){
  post.samp <- generate(input.list$model.fit, data.frame(),
                        ~ c(input.list$link.functions$mu(th.mu),
                            input.list$link.functions$K(th.K),
                            input.list$link.functions$alpha(th.alpha),
                            input.list$link.functions$c_(th.c),
                            input.list$link.functions$p(th.p)), n.samples = n.samp)
  data.frame(mu = post.samp[1,],
             K = post.samp[2,],
             alpha = post.samp[3,],
             c = post.samp[4,],
             p = post.samp[5,])
}


#' Sample from the posterior of the ETAS parameters
#'
#' @param input.list structured input `list` with at least two elements:
#' \itemize{
#' \item `model.fit`: `bru` object used to sample the posterior of the ETAS parameters
#' \item `link.functions`: `list` of functions to convert the ETAS parameters from the INLA scale to the ETAS scale
#' }
#' @param n.samp The number of samples to draw from the posteriors for the plot
#' @param post.samp : `data.frame` with columns mu, K, alpha, c, p and rows corresponding to different posterior samples. When `NULL` the function samples the joint posterior distribution `n.samp` times. The default is `NULL`.
#' @return `list`: with elements
#' \itemize{
#' \item `post.samp.df`:`data.frame` of posterior samples with `nrow = n.samp` and columns `mu, K, alpha, c, p` corresponding to ETAS parameters. If `post.samp` is not `NULL` it returns `post.samp`}
#' \item `pair.plot`: `ggplot` object reporting the pair plot between parameters samples. It is obtained using the `ggpairs` function of the `Ggally` library
#' @export
#'
#' @examples
post_pairs_plot <- function(input.list, n.samp, post.samp = NULL){
  if(is.null(input.list$model.fit)){
    stop('model.fit is missing, please provide a fitted model as bru object')
  }
  if(is.null(post.samp)){
    post.samp <- generate(input.list$model.fit, data.frame(),
                          ~ c(input.list$link.functions$mu(th.mu),
                              input.list$link.functions$K(th.K),
                              input.list$link.functions$alpha(th.alpha),
                              input.list$link.functions$c_(th.c),
                              input.list$link.functions$p(th.p)), n.samples = n.samp)
    post.samp.df <- data.frame(mu = post.samp[1,],
                               K = post.samp[2,],
                               alpha = post.samp[3,],
                               c = post.samp[4,],
                               p = post.samp[5,])
  }
  pair.plot <- ggpairs(post.samp.df, labeller = label_parsed, lower=list(continuous='density'))
  return(list(post.samp.df = post.samp.df,
              pair.plot = pair.plot))
}


#######
## MN: DESCRIPTION: Calculate the number of events in a time interval T1 to T2 given imposed events and ETAS
## MN: Arguments:
## MN: Returns: number of events
#' Calculate the integral of the ETAS conditional intensity
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
#' @param link.functions `list` of functions to trasnform the parameters from the internal scale to the ETAS scale
#'
#' @return Integral of the ETAS conditional intensity between `T1` and `T2` with minimum magnitude `M0`.
#' @export
#'
#' @examples
lambda.N <- function(th.mu, th.K, th.alpha, th.c, th.p, T1, T2, M0, Ht,
                     link.functions){
  theta_etas <- list(mu = link.functions$mu(th.mu[1]),
                     K = link.functions$K(th.K[1]),
                     alpha = link.functions$alpha(th.alpha[1]),
                     c = link.functions$c_(th.c[1]),
                     p = link.functions$p(th.p[1]))


  return( theta_etas$mu*(T2 - T1) + sum(exp(log.Lambda.h(theta = theta_etas,
                                                         th = Ht$ts,
                                                         mh = Ht$magnitudes,
                                                         M0 = M0,
                                                         T1 = T1, T2 = T2))) )
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
#'
#' @return A `list` of three objects:
#' \itemize{
#' \item `post.df`: `data.frame` containing posterior informations on the posterior distribution of the number of events
#' \item `post.plot` : `ggplot` object showing the posterior distribution of the expected number of events
#' \item `post.plot.shaded` : `ggplot` object showing the posterior distribution of the expected number of events, shaded region corresponds to the 0.025 and 0.975 quantiles of the distribution of the distribution of the number of events.}
#' @export
#'
#' @examples
get_posterior_N <- function(input.list){
  lambda.N.post <- predict(input.list$model.fit, # model fit
                           data.frame(), # data (empty because the history of the process is passed to the function below directly)
                           ~ lambda.N(th.mu, th.K, th.alpha, th.c, th.p,
                                      input.list$T12[1], input.list$T12[2], input.list$M0,
                                      input.list$catalog.bru,
                                      input.list$link.functions)) # target function
  N.seq <- floor(lambda.N.post$q0.025 - lambda.N.post$q0.025*0.10):ceiling(lambda.N.post$q0.975 + lambda.N.post$q0.975*0.10)
  N.post.df <- predict(input.list$model.fit, data.frame(),
                       ~ data.frame(N = N.seq,
                                    pdf = dpois(N.seq,
                                                lambda.N(th.mu, th.K,
                                                         th.alpha, th.c, th.p,
                                                         input.list$T12[1], input.list$T12[2], input.list$M0,
                                                         input.list$catalog.bru,
                                                         input.list$link.functions))))
  N.post.plot <- ggplot(N.post.df, aes(x = N, y = mean)) +
    geom_line() +
    geom_vline(xintercept = nrow(input.list$catalog.bru), linetype = 2) +
    ylab('pdf')

  N.post.plot.shaded <- N.post.plot +
    geom_ribbon(aes(xmax = N, xmin = N, ymin = q0.025, ymax = q0.975), alpha = 0.2,
                fill = 'blue')
  list(post.df = N.post.df,
       post.plot = N.post.plot,
       post.plot.shaded = N.post.plot.shaded)
}
