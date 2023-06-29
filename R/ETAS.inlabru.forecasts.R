#' Title
#'
#' @param post.samp a `data.frame` containing samples from the posterior distribution of ETAS parameters.
#' Each row of the `data.frame` corresponds to a different sample and the parameters are in the order
#' \eqn{\mu}{mu}, \eqn{K}{K}, \eqn{\alpha}{alpha}, \eqn{c}{c}, \eqn{p}{p}.
#' @param n.cat number of synthetic catalogues composing the forecast. If `n.cat` is greater than `nrow(post.samp)`, then,
#' `n.cat` rows are sampled uniformly and with replacement from the `post.samp`. If `n.cat` is smaller than `nrow(post.samp)`, then,
#' `n.cat` rows are sampled uniformly and without replacement from the `post.samp`. If `n.cat` is `NULL` or equal to `nrow(post.samp)`, `post.samp` is used at it is and `nrow(post.samp)` catalogues are generated.
#' @param beta.p parameter of the magnitude distribution
#' @param M0 cutoff magnitude, all the synthetic events will have magnitude greater than this value.
#' @param T1 starting time of the forecast
#' @param T2 end time of the forecast
#' @param Ht set of known events
#' @param ncore number of cores to be used to generate the synthetic catalogues in parallel.
#'
#' @return a `list` with two elements: `fore.df` is a `data.frame` containing all the synthetic catalogues composing the forecast.
#' The `data.frame` has four columns, `ts` for the occurrence time, `magnitudes` for the magnitude, `gen` with the generation of the event, and `cat.idx` with the catalogue identifier
#' The second element of the output `list` is `n.cat` which is the number of synthetic catalogues generated.
#' @seealso [generate_temporal_ETAS_synthetic()]
#' @export
Temporal.ETAS.forecast <- function(post.samp, n.cat, beta.p, M0, T1, T2, Ht, ncore = 1){
  if(n.cat > nrow(post.samp)){
    post.samp <- post.samp[sample(seq_len(nrow(post.samp)), n.cat, replace = TRUE),]
  } else if(n.cat < nrow(post.samp)){
    post.samp <- post.samp[sample(seq_len(nrow(post.samp)), n.cat),]
  }
  n.cat <- nrow(post.samp)
  synth.cat.list <- parallel::mclapply(seq_len(n.cat), \(x)
                                       generate_temporal_ETAS_synthetic(theta = post.samp[x,],
                                                                        beta.p = beta.p,
                                                                        M0 = M0,
                                                                        T1 = T1,
                                                                        T2 = T2,
                                                                        Ht = Ht),
                                       mc.cores = ncore)
  synth.cat.list.df <- lapply(synth.cat.list, \(x) do.call(rbind, x))
  synth.cat.with.events <- vapply(synth.cat.list.df, \(x) nrow(x) > 0, TRUE)
  idx.cat.with.events <- seq_len(n.cat)[synth.cat.with.events]
  if(length(idx.cat.with.events) == 0){
    df.out <- data.frame(ts = 0, magnitudes = 0, gen = 0, cat.idx = 0)
    return(list(fore.df = df.out[-1,],
                n.cat = n.cat))
  } else {
    # set catalogue identifier
    synth.cat.list.df <- lapply(idx.cat.with.events, \(x) cbind(synth.cat.list.df[[x]],
                                                                cat.idx = x))
    return(list(fore.df = do.call(rbind, synth.cat.list.df),
                n.cat = n.cat))
  }

}
