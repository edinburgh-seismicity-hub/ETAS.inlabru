
#' Generate summary information on the fitted ETAS model
#'
#' @param input.list Which has combined the input file (for link functions) and bru output (for marginals)
#'
#' @return Data frame summary and summary plot
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


#' Function to return a many (n.samp) samples from the posterior of the parameters
#'
#' @param input.list Which has combined the input file (for link functions) and bru output (for marginals)
#' @param n.samp The number of samples to draw from the posteriors
#'
#' @return n.samp samples drawn from the posteriors.
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

#######
## MN: DESCRIPTION: Calculate the number of events in a time interval T1 to T2 given imposed events and ETAS
## MN: Arguments:
## MN: Returns: number of events
#' Title
#'
<<<<<<< HEAD
#' @param th.mu
#' @param th.K
#' @param th.alpha
#' @param th.c
#' @param th.p
#' @param T1
#' @param T2
#' @param M0
=======
#' @param th.mu Background rate, `mu`, on the internal parameter scale
#' @param th.K ETAS triggering parameter `K` on the internal parameter scale
#' @param th.alpha ETAS triggering parameter `alpha` on the internal parameter scale
#' @param th.c ETAS triggering parameter `c` on the internal parameter scale
#' @param th.p ETAS triggering parameter `p` on the internal parameter scale
#' @param T1 Start of temporal model domain.
#' @param T2 End of temporal model domain.
#' @param M0 Minimum magnitude threshold
>>>>>>> 4bbeeb6f841c44ff042a96c001abbec873906594
#' @param Ht
#' @param link.functions
#'
#' @return
#' @export
#'
#' @examples
lambda.N <- function(th.mu, th.K, th.alpha, th.c, th.p, T1, T2, M0, Ht,
                     link.functions){
  theta_etas <- c(link.functions$mu(th.mu[1]),
                  link.functions$K(th.K[1]),
                  link.functions$alpha(th.alpha[1]),
                  link.functions$c_(th.c[1]),
                  link.functions$p(th.p[1]))


  return( theta_etas[1]*(T2 - T1) + sum(exp(log.Lambda_h2(th = theta_etas,
                                                          ti = Ht$ts,
                                                          mi = Ht$magnitudes,
                                                          M0 = M0,
                                                          T1 = T1, T2 = T2))) )
}

#######
## MN: DESCRIPTION: Samples the distribution of distributions of N
## MN: Arguments: input.list combining input file with bru model output
## MN: Returns: (i) DataFrame with summary statistics for N
##              (ii) Plot of the expected number of events
##              (iii) Plot with shaded stack of all distributions of N and a normal about the mean model

#' Title
#'
<<<<<<< HEAD
#' @param input.list
=======
#' @param input.list Which has combined the input file (for link functions) and bru output (for marginals)
>>>>>>> 4bbeeb6f841c44ff042a96c001abbec873906594
#'
#' @return
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
