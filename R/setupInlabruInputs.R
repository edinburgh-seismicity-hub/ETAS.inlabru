#' @title Create input list for ETAS Hawkes temporal model with catalogue
#' @decription
#' Function to create a default input file for the ETAS Hawkes temporal model
#' where a catalogue is specified in the input file.
#'
#' @param input_path path of the `txt` file containing experiment's information
#' @param num.threads Optional argument for the number of threads to be used by
#' parallel processing in inlabru/INLA
#'
#' @return The formatted input.list with the elements required for the temporal
#'   Hawkes model
#' @export
create_input_list_temporal_withCatalogue <- function(input_path,
                                                     num.threads = NULL) {
  con <- file(input_path)
  on.exit(close(con))
  par.raw <- readLines(con)
  for (i in seq_len(length(par.raw))) {
    row.i <- par.raw[[i]]
    if (grepl("start.date", row.i)) {
      # Do explicit assignment of the eval-result, to avoid package check
      # warnings
      start.date <- eval(parse(text = row.i))
    } else if (grepl("end.date", row.i)) {
      end.date <- eval(parse(text = row.i))
    } else if (grepl("magnitude.completeness", row.i)) {
      magnitude.completeness <- eval(parse(text = row.i))
    } else if (grepl("min.longitude", row.i)) {
      min.longitude <- eval(parse(text = row.i))
    } else if (grepl("max.longitude", row.i)) {
      max.longitude <- eval(parse(text = row.i))
    } else if (grepl("min.latitude", row.i)) {
      min.latitude <- eval(parse(text = row.i))
    } else if (grepl("max.latitude", row.i)) {
      max.latitude <- eval(parse(text = row.i))
    } else if (grepl("catalog.path", row.i)) {
      catalog.path <- eval(parse(text = row.i))
    } else if (grepl("catalog.header", row.i)) {
      catalog.header <- eval(parse(text = row.i))
    } else if (grepl("catalog.sep", row.i)) {
      catalog.sep <- eval(parse(text = row.i))
    } else if (grepl("catalog.skip", row.i)) {
      catalog.skip <- eval(parse(text = row.i))
    } else if (grepl("catalog.colnames", row.i)) {
      catalog.colnames <- eval(parse(text = row.i))
    } else if (grepl("a_mu", row.i)) {
      a_mu <- eval(parse(text = row.i))
    } else if (grepl("b_mu", row.i)) {
      b_mu <- eval(parse(text = row.i))
    } else if (grepl("a_K", row.i)) {
      a_K <- eval(parse(text = row.i))
    } else if (grepl("b_K", row.i)) {
      b_K <- eval(parse(text = row.i))
    } else if (grepl("a_alpha", row.i)) {
      a_alpha <- eval(parse(text = row.i))
    } else if (grepl("b_alpha", row.i)) {
      b_alpha <- eval(parse(text = row.i))
    } else if (grepl("a_c", row.i)) {
      a_c <- eval(parse(text = row.i))
    } else if (grepl("b_c", row.i)) {
      b_c <- eval(parse(text = row.i))
    } else if (grepl("a_p", row.i)) {
      a_p <- eval(parse(text = row.i))
    } else if (grepl("b_p", row.i)) {
      b_p <- eval(parse(text = row.i))
    } else if (grepl("mu.init", row.i)) {
      mu.init <- eval(parse(text = row.i))
    } else if (grepl("K.init", row.i)) {
      K.init <- eval(parse(text = row.i))
    } else if (grepl("alpha.init", row.i)) {
      alpha.init <- eval(parse(text = row.i))
    } else if (grepl("c.init", row.i)) {
      c.init <- eval(parse(text = row.i))
    } else if (grepl("p.init", row.i)) {
      p.init <- eval(parse(text = row.i))
    } else if (grepl("max_iter", row.i)) {
      max_iter <- eval(parse(text = row.i))
    } else if (grepl("max_step", row.i)) {
      max_step <- eval(parse(text = row.i))
    } else if (grepl("coef.t", row.i)) {
      coef.t <- eval(parse(text = row.i))
    } else if (grepl("DELTA", row.i)) {
      DELTA <- eval(parse(text = row.i))
    } else if (grepl("Nmax", row.i)) {
      Nmax <- eval(parse(text = row.i))
    } else if (grepl("n.periods", row.i)) {
      n.periods <- eval(parse(text = row.i))
    } else if (grepl("period.length", row.i)) {
      period.length <- eval(parse(text = row.i))
    } else if (grepl("start.date.fore", row.i)) {
      start.date.fore <- eval(parse(text = row.i))
    } else if (grepl("magnitude.update", row.i)) {
      magnitude.update <- eval(parse(text = row.i))
    } else if (grepl("output.name", row.i)) {
      output.name <- eval(parse(text = row.i))
    } else if (grepl("injection_vol2activity", row.i)) {
      injection_vol2activity <- eval(parse(text = row.i))
    } else if (grepl("injection_activityDecayRate", row.i)) {
      injection_activityDecayRate <- eval(parse(text = row.i))
    } else if (grepl("scenario", row.i)) {
      scenario <- eval(parse(text = row.i))
    }
  }
  # loading catalog
  catalog <- read.table(catalog.path,
    header = catalog.header,
    sep = catalog.sep,
    skip = catalog.skip
  )
  if (!catalog.header) {
    colnames(catalog) <- catalog.colnames
  }
  missing_names <- c("time_string", "Lon", "Lat", "magnitudes")
  missing_names <- missing_names[!(missing_names %in% colnames(catalog))]
  if (length(missing_names) > 0) {
    missing_msg <- c(
      "time_string" = "observed times",
      "Lon" = "observed longitudes",
      "Lat" = "observed latitudes",
      "magnitudes" = "observed magnitudes"
    )[missing_names]
    stop(paste0(
      "Error(s) in the catalog column names:\n",
      paste0("Please set the column name of ", missing_msg, " equal to '",
        missing_names, "'",
        collapse = "\n"
      )
    ))
  }
  #
  # catalog preparation
  start.date <- as.POSIXct(
    start.date,
    tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS")
  )
  end.date <- as.POSIXct(
    end.date,
    tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS")
  )
  catalog <- catalog %>%
    dplyr::mutate(
      time_date = as.POSIXct(
        .data$time_string,
        tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS")
      )
    ) %>%
    dplyr::filter(
      .data$time_date >= start.date,
      .data$time_date <= end.date,
      .data$Lon >= min.longitude,
      .data$Lon <= max.longitude,
      .data$Lat >= min.latitude,
      .data$Lat <= max.latitude,
      .data$magnitudes >= magnitude.completeness
    ) %>%
    dplyr::mutate(
      time_diff = as.numeric(difftime(.data$time_date,
        start.date,
        units = "days"
      ))
    )
  cat("Finish loading & preparing catalog", "\n")

  # create data.frame for inlabru
  data.bru <- data.frame(
    ts = catalog$time_diff,
    magnitudes = catalog$magnitudes,
    idx.p = seq_len(nrow(catalog))
  )
  # set up time interval and magnitude of completeness
  T1 <- 0
  T2 <- ceiling(as.numeric(difftime(end.date, start.date, units = "days")))
  M0 <- magnitude.completeness

  # priors
  link.f <- list(
    mu = \(x) gamma_t(x, a_mu, b_mu),
    K = \(x) loggaus_t(x, a_K, b_K),
    alpha = \(x) unif_t(x, a_alpha, b_alpha),
    c_ = \(x) unif_t(x, a_c, b_c),
    p = \(x) unif_t(x, a_p, b_p)
  )

  # initial value - convert from ETAS scale to internal scale
  th.init <- list(
    th.mu = inv_gamma_t(mu.init, a_mu, b_mu),
    th.K = inv_loggaus_t(K.init, a_K, b_K),
    th.alpha = inv_unif_t(alpha.init, a_alpha, b_alpha),
    th.c = inv_unif_t(c.init, a_c, b_c),
    th.p = inv_unif_t(p.init, a_p, b_p)
  )

  # options for inlabru
  if (is.null(max_step)) {
    bru.opt.list <- list(
      bru_verbose = 3, # type of visual output
      bru_max_iter = max_iter, # maximum number of iterations
      num.threads = num.threads,
      # bru_method = list(max_step = 0.5),
      # inla.mode = "experimental", # type of inla algorithm
      bru_initial = th.init
    ) # parameters initial values
  } else {
    bru.opt.list <- list(
      bru_verbose = 3, # type of visual output
      bru_max_iter = max_iter, # maximum number of iterations
      bru_method = list(max_step = max_step),
      num.threads = num.threads,
      # inla.mode = "experimental", # type of inla algorithm
      bru_initial = th.init
    ) # parameters initial values
  }
  # output list
  if (scenario == "temporal") {
    list(
      catalog = catalog,
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
      # model.fit = fit_etas,
      n.periods = n.periods,
      period.length = period.length,
      start.date.fore = start.date.fore,
      magnitude.update = magnitude.update,
      output.name = output.name
    )
  } else if (scenario == "injection") {
    list(
      catalog = catalog,
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
      # model.fit = fit_etas,
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


#' @title Create input list for ETAS Hawkes temporal model without catalogue
#'
#' @decription
#' Function to create a default input list for the ETAS Hawkes temporal model
#' where no catalogue is specified in the input file
#'
#' @param input_path Input file and path as a string
#' @param num.threads Optional argument for the number of threads to be used by
#' parallel processing by inlabru/INLA
#'
#' @return The formatted input.list with the elements required for the temporal
#'   Hawkes model
#' @export
#'
#' @examples
#' create_input_list_temporal_noCatalogue(
#'   system.file(
#'     "extdata",
#'     "user_input_synthetic_noCatalogue.txt",
#'     package = "ETAS.inlabru"
#'   )
#' )
create_input_list_temporal_noCatalogue <- function(input_path,
                                                   num.threads = NULL) {
  con <- file(input_path)
  on.exit(close(con))
  par.raw <- readLines(con)

  for (i in seq_len(length(par.raw))) {
    row.i <- par.raw[[i]]

    if (grepl("a_mu", row.i)) {
      # Do explicit assignment of the eval-result, to avoid package check
      # warnings
      a_mu <- eval(parse(text = row.i))
    } else if (grepl("b_mu", row.i)) {
      b_mu <- eval(parse(text = row.i))
    } else if (grepl("a_K", row.i)) {
      a_K <- eval(parse(text = row.i))
    } else if (grepl("b_K", row.i)) {
      b_K <- eval(parse(text = row.i))
    } else if (grepl("a_alpha", row.i)) {
      a_alpha <- eval(parse(text = row.i))
    } else if (grepl("b_alpha", row.i)) {
      b_alpha <- eval(parse(text = row.i))
    } else if (grepl("a_c", row.i)) {
      a_c <- eval(parse(text = row.i))
    } else if (grepl("b_c", row.i)) {
      b_c <- eval(parse(text = row.i))
    } else if (grepl("a_p", row.i)) {
      a_p <- eval(parse(text = row.i))
    } else if (grepl("b_p", row.i)) {
      b_p <- eval(parse(text = row.i))
    } else if (grepl("mu.init", row.i)) {
      mu.init <- eval(parse(text = row.i))
    } else if (grepl("K.init", row.i)) {
      K.init <- eval(parse(text = row.i))
    } else if (grepl("alpha.init", row.i)) {
      alpha.init <- eval(parse(text = row.i))
    } else if (grepl("c.init", row.i)) {
      c.init <- eval(parse(text = row.i))
    } else if (grepl("p.init", row.i)) {
      p.init <- eval(parse(text = row.i))
    } else if (grepl("max_iter", row.i)) {
      max_iter <- eval(parse(text = row.i))
    } else if (grepl("max_step", row.i)) {
      max_step <- eval(parse(text = row.i))
    } else if (grepl("coef.t", row.i)) {
      coef.t <- eval(parse(text = row.i))
    } else if (grepl("DELTA", row.i)) {
      DELTA <- eval(parse(text = row.i))
    } else if (grepl("Nmax", row.i)) {
      Nmax <- eval(parse(text = row.i))
    } else if (grepl("n.periods", row.i)) {
      n.periods <- eval(parse(text = row.i))
    } else if (grepl("period.length", row.i)) {
      period.length <- eval(parse(text = row.i))
    } else if (grepl("magnitude.update", row.i)) {
      magnitude.update <- eval(parse(text = row.i))
    } else if (grepl("output.name", row.i)) {
      output.name <- eval(parse(text = row.i))
    } else if (grepl("scenario", row.i)) {
      scenario <- eval(parse(text = row.i))
    }
  }
  # priors
  link.f <- list(
    mu = \(x) gamma_t(x, a_mu, b_mu),
    K = \(x) loggaus_t(x, a_K, b_K),
    alpha = \(x) unif_t(x, a_alpha, b_alpha),
    c_ = \(x) unif_t(x, a_c, b_c),
    p = \(x) unif_t(x, a_p, b_p)
  )

  # initial value - convert from ETAS scale to internal scale
  th.init <- list(
    th.mu = inv_gamma_t(mu.init, a_mu, b_mu),
    th.K = inv_loggaus_t(K.init, a_K, b_K),
    th.alpha = inv_unif_t(alpha.init, a_alpha, b_alpha),
    th.c = inv_unif_t(c.init, a_c, b_c),
    th.p = inv_unif_t(p.init, a_p, b_p)
  )

  # options for inlabru
  if (is.null(max_step)) {
    bru.opt.list <- list(
      bru_verbose = 3, # type of visual output
      bru_max_iter = max_iter, # maximum number of iterations
      num.threads = num.threads,
      # bru_method = list(max_step = 0.5),
      # inla.mode = "experimental", # type of inla algorithm
      bru_initial = th.init
    ) # parameters initial values
  } else {
    bru.opt.list <- list(
      bru_verbose = 3, # type of visual output
      bru_max_iter = max_iter, # maximum number of iterations
      bru_method = list(max_step = max_step),
      num.threads = num.threads,
      # inla.mode = "experimental", # type of inla algorithm
      bru_initial = th.init
    ) # parameters initial values
  }

  return(
    list(
      catalog = NULL,
      catalog.bru = NULL,
      time.int = c(NULL, NULL),
      T12 = c("T1", " T2"),
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
      # model.fit = fit_etas,
      n.periods = n.periods,
      period.length = period.length,
      start.date.fore = NULL,
      magnitude.update = magnitude.update,
      output.name = output.name
    )
  )
}
