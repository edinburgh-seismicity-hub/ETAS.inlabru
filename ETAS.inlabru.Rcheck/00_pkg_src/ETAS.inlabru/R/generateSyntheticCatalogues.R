## DESCRIPTION: Generates a whole catalogue by using a set of known events
## Arguments
##          Ht: Ether (i) Historic events to precondition sequence or (ii) notable in sequence events (e.g. Cofiorito)
##                Dataframe[t_h, M_h]
## Returns: The catalogue Dataframe[t_i, M_i, gen_i] including the events in Ht that are within T1 to T2
##                gen_i = -1 for Ht
##                gen_i = 0 for daughters of Ht
##                gen_i = 1 for background events
##                gen_i = 2 for first order daughters of background or second order daughters historic
##                gen_i = increments for successive generations

#' Generates a sythetic catalogue using the ETAS model
#'
#' @param theta ETAS parameters `data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`.
#' @param beta.p  Slope of GR relation: beta = b ln(10).
#' @param M0  The minimum magnitude in the synthetic catalogue.
#' @param T1 The start time for the synthetic catalogue [days].
#' @param T2 The end time for the synthetic catalogue [days].
#' @param Ht A catalogue history to impose on the synthetic sequence.
#' @param ncore Integer number of compute cores to use.
#' @return A data.frame of the temporal catalogue with columns `[t_i, M_i, gen_i]`
#' where, `t_i` are the times,  `M_i` the magnitudes, `gen_i` includes information about the generation number
#' @examples
#' ## EXAMPLE 1: Generate a 1000 day synthetic ETAS catalogue
#'
#' generate.temporal.ETAS.synthetic( theta=data.frame(mu=0.1, K=0.089, alpha=2.29, c=0.11, p=1.08), beta.p=log(10), M0=2.5, T1=0, T2=1000 )
#'
#'
#' ## EXAMPLE 2: To generate a 1000 day catalogue including a M6.7 event on day 500
#'
#' Ht <- data.frame(ts=c(500), magnitudes=c(6.7))
#' generate.temporal.ETAS.synthetic( theta=data.frame(mu=0.1, K=0.089, alpha=2.29, c=0.11, p=1.08), beta.p=log(10), M0=2.5, T1=0, T2=1000, Ht=Ht )

generate.temporal.ETAS.synthetic <- function(theta, beta.p, M0, T1, T2,
                            Ht = NULL, ncore = 1){
  #sample.ETAS <- function(theta, beta.p, M0, T1, T2,

  # Check that the end time is greater than the start time
  if(T2 < T1){
    stop('Error - right-end of time interval greater than left-end')
  }

  #########
  ### Background events
  # Expected number of background events
  n.bkg <- rpois(1, theta$mu*(T2 - T1))

  #cat('Background : ', n.bkg, '\n')
  # if no background events are generated initialize an empty data.frame
  if(n.bkg == 0){
    bkg.df <- data.frame(ts = 1, magnitudes = 1, gen = 0)
    bkg.df <- bkg.df[-1,]
  } else {
    # sample bkg events
    # if no bk.field.list element is passed it assumes uniform background rate

    # otherwise it samples using the information provided

    bkg.df <- data.frame(ts = runif(n.bkg, T1, T2),
                         magnitudes = rexp(n.bkg, beta.p) + M0,
                         gen = 1)
  }

  #########
  ### Generate from imposed events listed in Ht and the background events

  ## MN: Combine imposed and background events and add 1st generation daughters of imposed event set
  # if known events are provided
  if(!is.null(Ht)){
    #### TODO : the function has to generate all the points.
    # sample a generation from the known events
    gen.from.past <- sample.temporal.ETAS.generation(theta, beta.p, Ht, M0, T1, T2, ncore)
    # if at least an aftershock is produced
    if(nrow(gen.from.past) > 0){
      # set generation
      gen.from.past$gen = 0
      # Merge first generation and background events
      Gen.list <- list(rbind(gen.from.past, bkg.df))
    }
    else{
      Gen.list <- list(bkg.df)
    }

  }
  else{
    Gen.list <- list(bkg.df)
  }
  # stop if we have no background events and no events generated from known observations
  if(nrow(Gen.list[[1]]) == 0){
    #print(exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
    #stop('No events generated - increase theta1')
    return(Gen.list)
  }

  ## MN: Generate the daughters of backgorund and 1st daughters of imposed event set
  # initialize flag and gen counter
  flag = TRUE
  gen = 1
  # this goes until the condition inside the loop is met
  while(flag){
    # set parents
    parents <- Gen.list[[gen]]
    #cat('Gen : ', nrow(parents), '\n')
    #print(c(T1,T2))
    #print(range(parents$ts))
    # generate aftershocks
    triggered <- sample.temporal.ETAS.generation(theta, beta.p, parents,
                                   M0, T1, T2, ncore)
    #print(nrow(triggered))
    # stop the loop if there are no more aftershocks
    if(nrow(triggered) == 0){
      flag = FALSE}
    else{
      # set generations
      triggered$gen = gen + 1
      # store new generation
      Gen.list[[gen + 1]] = triggered
      # update generation counter
      gen = gen + 1
    }
  }
  Ht.to.gen <- Ht[Ht$ts >= T1 & Ht$ts <= T2, ]
  Ht.to.gen$gen = -1
  return( append(list(Ht), Gen.list) )
}



#' Take all previous parent events from `Ht=data.frame[ts, magnitudes]` and generates their daughters events using the ETAS model
#'
#' @param theta ETAS parameters `data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`.
#' @param beta.p Slope of GR relation: beta = b ln(10).
#' @param Ht The set of parent events in the form `data.frame[ts, magnitudes]`
#' @param M0 The minimum earthquake magnitude in the synthetic catalogue.
#' @param T1 The start time for the synthetic catalogue [days].
#' @param T2 The end time for the synthetic catalogue [days].
#' @param ncore The number of compute cores to use
#'
#' @return Return one generation of daughters from the parents in `Ht` in the form `data.frame(t_i, M_i)`.
#' @export
#'
#' @examples
#' # The parents are specified in Ht
#' Ht <- data.frame(ts=c(500), magnitudes=c(6.7))
#' sample.temporal.ETAS.generation( theta=data.frame(mu=0.1, K=0.089, alpha=2.29, c=0.11, p=1.08), beta.p=log(10), M0=2.5, T1=0, T2=1000, Ht=Ht )
sample.temporal.ETAS.generation <- function(theta, beta.p, Ht, M0, T1, T2, ncore = 1){

  # number of parents
  n.parent <- nrow(Ht)
  # calculate the aftershock rate for each parent in history (i.e. mean number daughters)
  trig.rates <- exp(log.Lambda_h2(th = theta,
                                  ti = Ht$ts, mi = Ht$magnitudes,
                                  M0 = M0, T1 = T1, T2 = T2))
  # extract number of aftershock for each parent (Sample poisson to deterime no daughters this realisation)
  n.ev.v <- sapply(1:n.parent, function(x) rpois(1, trig.rates[x]))

  # if no aftershock has to be generated returns empty data.frame
  if(sum(n.ev.v) == 0){
    app <- data.frame(x = 1, y = 1, ts = 1, magnitudes = 1)
    app <- app[-1,]
    return(app)
  }

  # identify parent events with the number of aftershocks > 0
  idx.p <- which(n.ev.v > 0)

  #print(sample.triggered(theta.v, beta.p, Sigma, Chol.M, n.ev.v[idx.p[1]], Ht[idx.p[1],], T1, T2, M0, bdy, crsobj))
  # sample (in parallel) the aftershocks for each parent
  sample.list <- parallel::mclapply(idx.p, function(idx)
    sample.temoral.ETAS.daughters(theta = theta, beta.p = beta.p, th = Ht$ts[idx],
                     n.ev = n.ev.v[idx], M0, T1, T2), mc.cores = ncore)

  # bind the data.frame in the list and return
  sample.pts <- dplyr::bind_rows(sample.list)
  sample.pts
}


##
#' Generate a sample of new events `data.frame(t_i, M_i)` of length `n.ev` for one parent event occuring at time `t_h` using the ETAS model.
#'
#' @param theta ETAS parameters `data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`.
#' @param beta.p Slope of GR relation: beta = b ln(10).
#' @param th Time of parent event [days].
#' @param n.ev The number of events to be placed.
#' @param M0 Minimum magnitude in synthetic catalogue.
#' @param T1 Start time for synthetic catalogue [days].
#' @param T2 End time for synthetic catalogue [days].
#'
#' @return Generate a sample of new events `data.frame(t_i, M_i)` from one parent
#' @export
#'
#' @examples
sample.temoral.ETAS.daughters <- function(theta, beta.p, th, n.ev, M0, T1, T2){
  # if the number of events to be placed is zero returns an empty data.frame
  if(n.ev == 0){
    samp.points <- data.frame(x = 1, y = 1, ts = 1, magnitudes = 1)
    samp.points <- samp.points[-1,]
    return(samp.points)
  }
  else{

    # Generate the time sample
    samp.ts <- sample.temporal.ETAS.times(theta, n.ev, th, T2)

    # Generate the magnitude sample
    samp.mags <- sample.GR.magnitudes(n=n.ev, beta.p=beta.p, M0=M0)
    print(c(length(samp.mags), length(samp.ts)))

    # Combine to build output synthetic catalogue for single parent
    samp.points <- data.frame(ts = samp.ts, magnitudes = samp.mags)
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    return(samp.points[!is.na(samp.points$ts),])
  }
}


#' Return a sample of magnitudes drawn from the GR distribution
#'
#' @param n Number of events in the sample.
#' @param beta.p Related to the b-value via `b ln(10)`.
#' @param M0 Minimum magnitude for the sample.
#'
#' @return A list of magnitudes of length `n` drawn from a GR distribution.
#' @export
#'
#' @examples
#' sample.GR.magnitudes(n=100, beta.p=log(10), M0=2.5)
sample.GR.magnitudes <- function(n, beta.p, M0) {
  return(rexp(n, beta.p) + M0)
  }


#' Sampling times for events triggered by a parent at th according to the ETAS triggering function
#'
#' @param theta ETAS parameters `data.frame(mu=mu, K=K, alpha=alpha, c=c, p=p)`.
#' @param n.ev Number of events to return in the sample in time domain (th, T2].
#' @param th Time of the parent event producing n.ev daughters.
#' @param T2 End time of model domain.
#'
#' @return t.sample A list of times in the interval [0, T2] distributed according to the ETAS triggering function.
#' @export
#'
#' @examples
sample.temporal.ETAS.times <- function(theta, n.ev, th, T2){
  if(n.ev == 0){
    df <- data.frame(ts = 1, x = 1, y = 1, magnitudes = 1, gen = 0)
    return(df[-1,])
  }
  bound.l <- 0 #Int.ETAS.time.trig.function(th.p, th, T)
  bound.u <- Int.ETAS.time.trig.function(theta, th, T2)
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  t.sample <- Inv.Int.ETAS.time.trig.function(theta, unif.s, th)
  return( t.sample )
}


#############################
#### Injection rate function
## Forward time integrated function for exponential rate decay
#' Title
#'
#' @param a Event rate per unit volume injected
#' @param V.i Injected volume
#' @param tau Decau rate [days]
#' @param T.i Time of injection event
#' @param T2
#'
#' @return
#' @export
#'
#' @examples
IntInjecIntensity <- function(a=50, V.i=1, tau=10, T.i, T2){
  expected.injection.events <- - tau*V.i*a* ( exp(-(T2-T.i)/tau ) -1 )
  return(expected.injection.events)
}

## Returns end time given a ...
#' Title
#'
#' @param a Event rate per unit volume injected
#' @param V.i Injected volume
#' @param tau Decau rate [days]
#' @param T.i Time of injection event
#' @param number.injected.events
#'
#' @return
#' @export
#'
#' @examples
Inv.IntInjecIntensity <- function(a=50, V.i=1, tau=10, T.i, number.injected.events){
  endTime <- T.i - tau*log(1 - number.injected.events / (tau*V.i*a))
  return(endTime)
}

#' Title
#'
#' @param a Induced event rate per unit volume.
#' @param V.i Injected volume
#' @param tau Decay rate [days].
#' @param beta.p Related to the b-value via `b ln(10)`.
#' @param M0 Minimum magnitude threshold.
#' @param T.i Time of injection [days].
#' @param T2 End of temporal model domain [days].
#'
#' @return Catalogue of parent events induced by injection data.frame(times, magnitudes)
#' @export
#'
#' @examples
sample.temoral.injection.events <- function(a=50, V.i=1, tau=10, beta.p, M0, T.i, T2){
  bound.l <- 0 #It(th.p, th, T)
  bound.u <- IntInjecIntensity(a=a, V.i=V.i, tau=tau, T.i=T.i, T2=T2)
  n.ev <- rpois( 1, bound.u  )
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  sample.ts <- Inv.IntInjecIntensity(a=a, V.i=V.i, tau=tau, T.i=T.i, number.injected.events=unif.s)

  samp.mags <- rexp(n.ev, rate = beta.p) + M0

  samp.points <- data.frame(ts = sample.ts, magnitudes = samp.mags)
  return(samp.points[!is.na(samp.points$ts),])
}

