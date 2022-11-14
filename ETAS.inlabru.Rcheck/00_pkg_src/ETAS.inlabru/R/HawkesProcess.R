## MN: DESCRIPTION: Fit a model according to the information in the input list
## MN: Arguments: input list containing parameters
## MN: Returns: bru output for the fitted ETAS model as a dataframe
#' Fits the remporal ETAS model and returns the results. This function decomposes the input.list for the `Hawkes.bru2`` function.
#'
#' @param input.list All input data and parameters are passed to inlabru via this structured list.
#'
#' @return The inlabru Hawkes process results data.frame
#' @export
#'
#' @examples
Temporal.ETAS.fit <- function(input.list){
  cat('Start model fitting', '\n')
  fit_etas <- Hawkes.bru2(sample.s = input.list$catalog.bru, # data
                          M0 = input.list$M0, # magnitude of completeness
                          T1 = input.list$T12[1],
                          T2 = input.list$T12[2], # time domain
                          link.functions = input.list$link.functions, # link functions
                          coef.t. = input.list$coef.t, # binning parameter (delta)
                          delta.t. = input.list$delta.t, # binning parameter (Delta)
                          N.max. = input.list$Nmax, # binning parameter (n.max)
                          bru.opt = input.list$bru.opt.list) # bru options
  cat('Finish model fitting', '\n')
  fit_etas
}


#######
## MN: DESCRIPTION: Function to fit an ETAS Hawkes process model to catalogue data
## MN: Arguments:
## MN: Returns: The inlabru Hawkes process results data.frame
#' Function to fit Hawkes process model
#'
#' @param sample.s
#' @param M0 Minimum magnitude threshold
#' @param T1 Start of temporal model domain [days].
#' @param T2 End of temporal model domain [days].
#' @param link.functions
#' @param coef.t. TimeBinning parameter:
#' @param delta.t. TimeBinning parameter:
#' @param N.max. TimeBinning parameter: Number of bins
#' @param bru.opt Runtime options for inlabru: See https://inlabru-org.github.io/inlabru/reference/bru_call_options.html
#'
#' @return The inlabru Hawkes process results data.frame
#' @export
#'
#' @examples
Hawkes.bru2 <- function(sample.s, M0, T1, T2, link.functions = NULL,
                        coef.t., delta.t., N.max., bru.opt){

  # Expected number of background events
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  # this is the expression of log(Lambda0)

  ## Create the grid for the XX integration
  cat('Start creating grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    time.grid(data.point = sample.s[idx,],
              coef.t = coef.t.,
              delta.t = delta.t.,
              T2. = T2, N.exp. = N.max.
    )
  }
  df.j$counts <- 0
  df.j$exposures <- 1
  df.j$part = 'triggered'

  t.names <- unique(df.j$t.ref_layer)
  time.sel <- df.j[vapply(t.names, \(bname) match(TRUE, df.j$t.ref_layer == bname), 0L), , drop = FALSE]
  Imapping <- match(df.j$t.ref_layer, t.names)

  cat('Finished creating grid, time ', Sys.time() - time.g.st, '\n')

  # MN: Define local function to calculate log-likelihood triggered contribution of one event to each bin.
  # MN: h is to denote it is from past events
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, list.input_, ncore_ = ncore){
    theta_ <- c(0,
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))

    #cat('theta - LogL', theta_, '\n')
    comp. <- compute.grid(param. = theta_, list.input_ = list.input_)
    #print(sum(is.na(comp.list$It)))
    #print(sum(is.infinite(comp.list$It)))
    out <- theta_[3]*(list.input_$df_grid$magnitudes - list.input_$M0) + log(theta_[2] + 1e-100) + log(comp. + 1e-100)
    out
  }

  # creating formula for past events contributions to integrated lambda
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')

  ## MN: Function to calculate sum of log intensities using past events
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, M0){

    ## MN: Rescale parameters to internal parameter scale
    if(is.null(link.functions)){
      th.p <- c(th.mu[1], th.K[1], th.alpha[1], th.c[1], th.p[1])
    }
    else{
      th.p <- c(link.functions$mu(th.mu[1]),
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))
    }

    ## MN: Parallel looping over the historic event magnitudes and times QQ why mean?  ALSO past in grid or historic?
    ## MN: Finn's trick for making more stable - QQ is this to avoid large numbers??
    mean(unlist(mclapply(tt, \(x) {
      th_x <- th < x
      log(lambda_2(th = th.p, t = x, ti.v = th[th_x],
                   mi.v = mh[th_x], M0 = M0))
    },
    mc.cores = 5)))
  }

  list.input <- list(df_grid = df.j, M0 = M0, Imapping = Imapping, time.sel = time.sel,
                     sample.s = sample.s)
  data.input = bind_rows(df.0, df.s, df.j)      ## Combine
  list.input <- append(list.input,
                       list(idx.bkg = data.input$part == 'background',
                            idx.trig = data.input$part == 'triggered',
                            idx.sl = data.input$part == 'SL'))
  predictor.fun <- function(th.mu, th.K, th.alpha, th.c, th.p,
                            list.input, T1, T2, M0){

    out <- rep(0, nrow(data.input))
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) + log(T2 - T1)
    out[list.input$idx.trig] <- logLambda.h.inla(th.K = th.K, th.alpha = th.alpha,
                                                 th.c = th.c, th.p = th.p,
                                                 list.input_ = list.input)
    out[list.input$idx.sl] <- loglambda.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                             th.c = th.c, th.p = th.p,
                                             tt = list.input$sample.s$ts,
                                             th = list.input$sample.s$ts,
                                             mh = list.input$sample.s$magnitudes,
                                             M0 = M0)
    out
  }

  merged.form <- counts ~ predictor.fun(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                        th.c = th.c, th.p = th.p, list.input = list.input,
                                        T1= T1, T2 = T2, M0 = M0)
  cmp.part <- counts ~ -1 +
    th.mu(1, model = 'linear', mean.linear = 0 , prec.linear = 1) +
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1)

  return( bru(formula = merged.form, components = cmp.part, data = data.input, family = 'Poisson',
      options = append(bru.opt, list(E = data.input$exposure))) )

}


