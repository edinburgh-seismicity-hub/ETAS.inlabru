input.file.to.list <- function(input_path){
  con <- file(input_path)
  on.exit(close(con))
  par.raw <- readLines(con)
  for(i in 1:length(par.raw)){
    row.i <- par.raw[[i]]
    if(grepl('start.date', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('end.date', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.completeness', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.longitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max.longitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.latitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max.latitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.path', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.header', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.sep', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.skip', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.colnames', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('mu.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('K.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('alpha.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('c.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('p.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_iter', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_step', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('coef.t', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('DELTA', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('Nmax', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('n.periods', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('period.length', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('start.date.fore', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.update', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('output.name', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('injection_vol2activity', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('injection_activityDecayRate', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('scenario', row.i)){
      eval(parse(text = row.i))
    }
  }
  # loading catalog
  catalog <- read.table(catalog.path,
                        header = catalog.header,
                        sep = catalog.sep,
                        skip = catalog.skip)
  if(!catalog.header){
    colnames(catalog) <- catalog.colnames
  }
  if(!('time_string' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed time equal to "time_string"')
  }
  if(!('Lon' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed longitudes equal to "time_string"')
  }
  if(!('Lat' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed latitude equal to "time_string"')
  }
  if(!('magnitudes' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed magnitudes equal to "time_string"')
  }
  #
  # catalog preparation
  start.date <- as.POSIXct(start.date, format = '%Y-%m-%d %H:%M:%OS')
  end.date <- as.POSIXct(end.date)
  catalog <- catalog %>%
    mutate(time_date = as.POSIXct(gsub(pattern = 'T',
                                       replacement = ' ',
                                       x = time_string),
                                  format = '%Y-%m-%d %H:%M:%OS')) %>%
    filter(time_date >= start.date,
           time_date <= end.date,
           Lon >= min.longitude,
           Lon <= max.longitude,
           Lat >= min.latitude,
           Lat <= max.latitude,
           magnitudes >= magnitude.completeness) %>%
    mutate(time_diff = as.numeric(difftime(time_date, start.date, units = 'days')))
  cat('Finish loading & preparing catalog', '\n')

  # create data.frame for inlabru
  data.bru <- data.frame(ts = catalog$time_diff,
                         magnitudes = catalog$magnitudes,
                         idx.p = seq_len(nrow(catalog)))
  # set up time interval and magnitude of completeness
  T1 <- 0
  T2 <- ceiling(as.numeric(difftime(end.date, start.date, units = 'days')))
  M0 <- magnitude.completeness

  # priors
  link.f <- list(mu = \(x) gamma.t(x, a_mu, b_mu),
                 K = \(x) loggaus.t(x, a_K, b_K),
                 alpha = \(x) unif.t(x, a_alpha, b_alpha),
                 c_ = \(x) unif.t(x, a_c, b_c),
                 p = \(x) unif.t(x, a_p, b_p))

  # initial value - convert from ETAS scale to internal scale
  th.init <- list(th.mu = inv.gamma.t(mu.init, a_mu, b_mu),
                  th.K = inv.loggaus.t(K.init, a_K, b_K),
                  th.alpha = inv.unif.t(alpha.init, a_alpha, b_alpha),
                  th.c = inv.unif.t(c.init, a_c, b_c),
                  th.p = inv.unif.t(p.init, a_p, b_p) )

  # options for inlabru
  if(is.null(max_step)){
    bru.opt.list <- list(bru_verbose = 3, # type of visual output
                         bru_max_iter = max_iter, # maximum number of iterations
                         num.threads = 5,
                         #bru_method = list(max_step = 0.5),
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
  } else {
    bru.opt.list <- list(bru_verbose = 3, # type of visual output
                         bru_max_iter = max_iter, # maximum number of iterations
                         bru_method = list(max_step = max_step),
                         num.threads = 5,
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values

  }
  # output list
  if(scenario=="temporal"){
    list(catalog = catalog,
         catalog.bru = data.bru,
         time.int = c(start.date, end.date),
         T12 = c(T1, T2),
         lat.int = c(min.latitude, max.latitude),
         lon.int = c(min.longitude, max.longitude),
         M0 = M0,
         link.functions = link.f,
         bru.opt.list = bru.opt.list,
         coef.t = coef.t,
         delta.t = DELTA,
         Nmax = Nmax,
         #model.fit = fit_etas,
         n.periods = n.periods,
         period.length = period.length,
         start.date.fore = start.date.fore,
         magnitude.update = magnitude.update,
         output.name = output.name
    )
  } else if(scenario=="injection") {
    list(catalog = catalog,
         catalog.bru = data.bru,
         time.int = c(start.date, end.date),
         T12 = c(T1, T2),
         lat.int = c(min.latitude, max.latitude),
         lon.int = c(min.longitude, max.longitude),
         M0 = M0,
         link.functions = link.f,
         bru.opt.list = bru.opt.list,
         coef.t = coef.t,
         delta.t = DELTA,
         Nmax = Nmax,
         #model.fit = fit_etas,
         n.periods = n.periods,
         period.length = period.length,
         start.date.fore = start.date.fore,
         magnitude.update = magnitude.update,
         output.name = output.name,
         injection_vol2activity = injection_vol2activity,
         injection_activityDecayRate = injection_activityDecayRate
    )
  }
}


createDeafultInputList.temporal.noCatalogue <- function(input_path){

  con <- file(input_path)
  on.exit(close(con))
  par.raw <- readLines(con)

  for(i in 1:length(par.raw)){
    row.i <- par.raw[[i]]

    if(grepl('a_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('mu.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('K.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('alpha.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('c.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('p.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_iter', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_step', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('coef.t', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('DELTA', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('Nmax', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('n.periods', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('period.length', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.update', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('output.name', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('scenario', row.i)){
      eval(parse(text = row.i))
    }
  }
  # priors
  link.f <- list(mu = \(x) gamma.t(x, a_mu, b_mu),
                 K = \(x) loggaus.t(x, a_K, b_K),
                 alpha = \(x) unif.t(x, a_alpha, b_alpha),
                 c_ = \(x) unif.t(x, a_c, b_c),
                 p = \(x) unif.t(x, a_p, b_p))

  # initial value - convert from ETAS scale to internal scale
  th.init <- list(th.mu = inv.gamma.t(mu.init, a_mu, b_mu),
                  th.K = inv.loggaus.t(K.init, a_K, b_K),
                  th.alpha = inv.unif.t(alpha.init, a_alpha, b_alpha),
                  th.c = inv.unif.t(c.init, a_c, b_c),
                  th.p = inv.unif.t(p.init, a_p, b_p) )

  # options for inlabru
  if(is.null(max_step)){
    bru.opt.list <- list(bru_verbose = 3, # type of visual output
                         bru_max_iter = max_iter, # maximum number of iterations
                         num.threads = 5,
                         #bru_method = list(max_step = 0.5),
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
  } else {
    bru.opt.list <- list(bru_verbose = 3, # type of visual output
                         bru_max_iter = max_iter, # maximum number of iterations
                         bru_method = list(max_step = max_step),
                         num.threads = 5,
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
  }

  return(
    list(catalog = NULL,
         catalog.bru = NULL,
         time.int = c(NULL,NULL),
         T12 = c('T1',' T2'),
         lat.int = c(-90, 90),
         lon.int = c(-180, 180),
         M0 = NULL,
         mu.init = mu.init,
         K.init = K.init,
         alpha.init = alpha.init,
         c.init = c.init,
         p.init = p.init,
         a_mu = a_mu,
         b_mu = b_mu,
         a_K = a_K,
         b_K = b_K,
         a_alpha = a_alpha,
         b_alpha = b_alpha,
         a_c = a_c,
         b_c = b_c,
         a_p = a_p,
         b_p = b_p,
         max_iter = max_iter,
         max_step = max_step,
         link.functions = link.f,
         bru.opt.list = bru.opt.list,
         coef.t = coef.t,
         delta.t = DELTA,
         Nmax = Nmax,
         #model.fit = fit_etas,
         n.periods = n.periods,
         period.length = period.length,
         start.date.fore = NULL,
         magnitude.update = magnitude.update,
         output.name = output.name
    )
  )
}
