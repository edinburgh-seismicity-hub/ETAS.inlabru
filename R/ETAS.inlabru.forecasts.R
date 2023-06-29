#' Title
#'
#' @param post.samp a `data.frame` containing samples from the posterior distribution of ETAS parameters.
#' Each row of the `data.frame` corresponds to a different sample and the parameters are in the order $\mu, K, \alpha, c, p$.
#' @param n.cat number of synthetic catalogues composing the forecast. If `n.cat` is greater than `nrow(post.samp)`, then,
#' `n.cat` rows are sampled uniformly and with replacement from the `post.samp`. If `n.cat` is smaller than `nrow(post.samp)`, then,
#' `n.cat` rows are sampled uniformly and without replacement from the `post.samp`. If `n.cat` is `NULL` or equal to `nrow(post.samp)`, `post.samp` is used at it is and `nrow(post.samp)` catalogues are generated.
#' @param beta.p parameter of the magnitude distribution
#' @param M0 cutoff magnitude, all the synthetic events will have magnitude greater than this value.
#' @param T1 starting time of the forecast
#' @param T2 end time of the forecast
#' @param Ht set of known events
#' @param ncore number of cores to be used to generate the syntehtic catalogues in parallel.
#'
#' @return a `data.frame` containing all the synthetic catalogues composing the forecast.
#' The `data.frame` has four columns: `ts` for the occurence time, `magnitudes` for the magnitude, `gen` with the generation of the event, and `cat.idx` with the catalogue identifier
#' @seealso [generate.temporal.ETAS.synthetic()]
#' @export
Temporal.ETAS.forecast <- function(post.samp, n.cat, beta.p, M0, T1, T2, Ht, ncore = 1){
  if(n.cat > nrow(post.samp)){
    post.samp <- post.samp[sample(seq_len(nrow(post.samp)), n.cat, replace = TRUE),]
  } else if(n.cat < nrow(post.samp)){
    post.samp <- post.samp[sample(seq_len(nrow(post.samp)), n.cat),]
  }
  synth.cat.list <- parallel::mclapply(seq_len(nrow(post.samp)), \(x)
                                       generate.temporal.ETAS.synthetic(theta = post.samp[x,],
                                                                        beta.p = beta.p,
                                                                        M0 = M0,
                                                                        T1 = T1,
                                                                        T2 = T2,
                                                                        Ht = Ht),
                                       mc.cores = ncore)
  synth.cat.list.df <- lapply(synth.cat.list, \(x) do.call(rbind, x))
  # set catalogue identifier
  synth.cat.list.df <- lapply(seq_len(n.cat), \(x) cbind(synth.cat.list.df[[x]],
                                                         cat.idx = x))
  do.call(rbind, synth.cat.list.df)
}
