
library(INLA)
library(inlabru)
library(hypergeo)
library(stats)
library(SuppDists)
library(foreach)

########## 1- defining new auxiliary functions ##########

# finding mainshock
mainshock <- function(list.input){
  subset(list.input$catalog.bru , magnitudes==max(list.input$catalog.bru$magnitudes))
}


# calculating temporal Mc and adding it to input.list
assign_Mc <- function(list.input, M0, G, H, b) {
  for (i in 1:nrow(list.input$catalog.bru)){
    if (list.input$catalog.bru$ts[i] < mainshock(list.input)$ts){
      list.input$catalog.bru$Mc[i] <- M0
    }
    else if (list.input$catalog.bru$ts[i] > mainshock(list.input)$ts){
      Mc_t <-  mainshock(list.input)$magnitudes - G - H * log10(list.input$catalog.bru$ts[i]-mainshock(list.input)$ts)
      list.input$catalog.bru$Mc[i] <- max(Mc_t, M0)
    }
  }
  return(list.input)
}



########## 2- modifying functions for linearisation ##########

# ETAS triggering function
modified_gt <- function(theta, t, th, mh, M0, G, H, b){
  if(is.list(theta)){
    mu <- theta$mu
    alpha <- theta$alpha
    K <- theta$K
    c <- theta$c
    p <- theta$p
  } else {
    mu <- theta[1]
    K <- theta[2]
    alpha <- theta[3]
    c <- theta[4]
    p <- theta[5]
  }

  output <- rep(0,length(th))
  t.diff <- t - th
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    Mainshock <- mainshock(list.input)
    Mm <- Mainshock$magnitudes
    Tm <- Mainshock$ts
    incompletness_end_time <- Tm + 10^((Mm-M0-G)/H)
    mct <- ifelse(t>Tm & t <= incompletness_end_time, Mm-G-(H*log10(t- Tm)), M0)
    log.out <- log(K) + alpha*(mh-M0) - p*log(1 + t.diff[!neg]/c) + log(10^(-b*(mct-M0)))
    output[!neg] <- exp(log.out)
  }
  else{
    output
  }
  output
}


# ETAS conditional intensity
modified_cond.lambda <- function(theta, t, th, mh, M0, G, H, b){
  if(is.null(th) | all(th > t)){
    theta$mu
  }
  theta$mu + sum(modified_gt(theta, t, th, mh, M0, G, H, b))
}


# modified_log.Lambda.h <- function(){}  # Not necessary now!



########## 3- modifying functions for binning ##########

# Computing the integral of the Omori's law
modified_It_df <- function(param_, time.df, M0, G, H, b){
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_c <- param_[4]
  param_p <- param_[5]
  T.l <- pmax(tth, T1b)

  Mainshock <- mainshock(list.input)
  Mm <- Mainshock$magnitudes
  Tm <- Mainshock$ts
  incompletness_end_time <- Tm + 10^((Mm-M0-G)/H)

  fun.l <- (1 + (T.l-tth)/param_c)^(1-param_p)
  fun.u <- (1 + (T2b-tth)/param_c)^(1-param_p)
  scenario_1 <- (param_c/(param_p-1)) * (fun.l - fun.u)

  fun.uu <- (param_c/(1-param_p)) * (((T2b-tth)/param_c+1)^(1-param_p))  * (10^(b*(G+M0-Mm))) * ((T2b-Tm)^(b*H))  * (((Tm-T2b)/(Tm-tth+param_c))^(-b*H))  * hypergeo(-b*H , 1-param_p , 2-param_p , -(T2b-tth+param_c)/(tth-Tm-param_c))
  fun.ll <- (param_c/(1-param_p)) * (((T.l-tth)/param_c+1)^(1-param_p))  * (10^(b*(G+M0-Mm))) * ((T.l-Tm)^(b*H))  * (((Tm-T.l)/(Tm-tth+param_c))^(-b*H))  * hypergeo(-b*H , 1-param_p , 2-param_p , -(T.l-tth+param_c)/(tth-Tm-param_c))
  scenario_2 <- Re(fun.uu - fun.ll)

  #saveRDS(is.nan(scenario_2) , "FARNAZ_NaNs.rds")
  ifelse(tth>Tm & tth<=incompletness_end_time, scenario_2, scenario_1)
}


# Computing the integral of the Omori's law efficiently
modified_compute.grid <- function(param., list.input_, M0, G, H, b){
  It.vec <- modified_It_df(param_=param., time.df=list.input_$time.sel, M0=M0, G=G, H=H, b=b)
  It.vec[list.input_$Imapping]
}



########## 4- modifying internal functions in ETAS fitting ##########

# Temporal.ETAS.fit
modified_Temporal.ETAS.fit <- function(input.list){
  cat('Start model fitting', '\n')
  fit_etas <- modified_Temporal.ETAS(total.data = input.list$catalog.bru,
                                     M0 = input.list$M0,
                                     T1 = input.list$T12[1],
                                     T2 = input.list$T12[2],
                                     link.functions = input.list$link.functions,
                                     coef.t. = input.list$coef.t,
                                     delta.t. = input.list$delta.t,
                                     N.max. = input.list$Nmax,
                                     bru.opt = input.list$bru.opt.list,
                                     G = input.list$G,
                                     H = input.list$H,
                                     b = input.list$b)
  cat('Finish model fitting', '\n')
  fit_etas
}



# Temporal.ETAS
# Things  to be changed = logLambda.h.inla / loglambda.inla / out / and some arguments
modified_Temporal.ETAS <- function(total.data, M0, T1, T2, link.functions = NULL,
                          coef.t., delta.t., N.max., bru.opt, G, H, b){

  if(sum(is.na(total.data))>0){
    print("Some input data is NA; removing rows")
    total.data <- na.omit(total.data)
  }

  idx.remove <- total.data$ts > T2
  if(sum(idx.remove) > 0){
    total.data <- total.data[total.data$ts > T2, ]
    warning('Removing events after T2')
  }
  idx.sample <- total.data$ts > T1 & total.data$ts < T2
  sample.s <- total.data[idx.sample, ]

  # Expected number of background events
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  # this is the expression of log(Lambda0)

  ## Create the grid for the XX integration
  cat('Start creating grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(total.data), .combine = rbind) %do% {
    time.grid(data.point = total.data[idx,],
              coef.t = coef.t.,
              delta.t = delta.t.,
              T1. = T1,
              T2. = T2,
              N.exp. = N.max.
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
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, list.input_, ncore_ = ncore, G, H, b){
    theta_ <- c(0,
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))

    #cat('theta - LogL', theta_, '\n')
    comp. <- modified_compute.grid(param. = theta_, list.input_ = list.input_, M0=M0, G=G, H=H, b=b)
    #print(sum(is.na(comp.list$It)))
    #print(sum(is.infinite(comp.list$It)))
    out <- theta_[3]*(list.input_$df_grid$magnitudes - list.input_$M0) + log(theta_[2] + 1e-100) + log(comp. + 1e-100)
    out
  }

  # creating formula for past events contributions to integrated lambda
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')

  ## MN: Function to calculate sum of log intensities using past events
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, M0, G, H, b){

    ## MN: Rescale parameters to internal parameter scale
    if(is.null(link.functions)){
      th.p <- list(mu = th.mu[1],
                   K = th.K[1], alpha = th.alpha[1],
                   c = th.c[1], p = th.p[1])
    }
    else{
      th.p <- list(mu = link.functions$mu(th.mu[1]),
                   K = link.functions$K(th.K[1]),
                   alpha = link.functions$alpha(th.alpha[1]),
                   c = link.functions$c_(th.c[1]),
                   p = link.functions$p(th.p[1]))
    }

    ## MN: Parallel looping over the historic event magnitudes and times QQ why mean?  ALSO past in grid or historic?
    ## MN: Finn's trick for making more stable - QQ is this to avoid large numbers??
    out <- mean(unlist(lapply(tt, \(x) {
      th_x <- th < x
      log(modified_cond.lambda(theta = th.p, t = x, th = th[th_x],
                      mh = mh[th_x], M0 = M0, G=G, H=H, b=b))
    })))#,
    #mc.cores = 5)))
    return(out)
  }

  list.input <- list(df_grid = df.j, M0 = M0, Imapping = Imapping, time.sel = time.sel,
                     sample.s = sample.s, total.data = total.data)
  data.input = bind_rows(df.0, df.s, df.j)      ## Combine
  list.input <- append(list.input,
                       list(idx.bkg = data.input$part == 'background',
                            idx.trig = data.input$part == 'triggered',
                            idx.sl = data.input$part == 'SL'))
  predictor.fun <- function(th.mu, th.K, th.alpha, th.c, th.p,
                            list.input, T1, T2, M0, G, H, b){

    out <- rep(0, nrow(data.input))
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) + log(T2 - T1)
    out[list.input$idx.trig] <- logLambda.h.inla(th.K = th.K, th.alpha = th.alpha,
                                                 th.c = th.c, th.p = th.p,
                                                 list.input_ = list.input,
                                                 G=G, H=H, b=b)
    out[list.input$idx.sl] <- loglambda.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                             th.c = th.c, th.p = th.p,
                                             tt = list.input$sample.s$ts,
                                             th = list.input$total.data$ts,
                                             mh = list.input$total.data$magnitudes,
                                             M0 = M0,
                                             G = G,
                                             H = H,
                                             b = b)
    out
  }

  merged.form <- counts ~ predictor.fun(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                        th.c = th.c, th.p = th.p, list.input = list.input,
                                        T1= T1, T2 = T2, M0 = M0, G = G, H = H, b = b)
  cmp.part <- counts ~ -1 +
    th.mu(1, model = 'linear', mean.linear = 0 , prec.linear = 1) +
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1)

  return( bru(formula = merged.form, components = cmp.part, data = data.input, family = 'Poisson',
              options = append(bru.opt, list(E = data.input$exposure))) )

}






