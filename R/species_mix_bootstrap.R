#' @title bootstrapping of species_mix or regional_mix fitted models.
#' @param object A species_mix model object
#' @param nboot The number of bootstraps to run.
#' @param type the type of bootstrap to use, options are "SimpleBoot" which is a parametric bootstrap,
#' or "BayesBoot" which uses a Bayesian Bootstrap. "BayesBoot" is the default
#' @param mc.cores The number of cores to use when running a bootstrap. The default is 1.
#' @param quiet If TRUE, do not print progress of bootstrap.
#' @param \\dots Ignored
#' @description This function can take a while to run -- it is a bootstrap function. nboot re-samples of the data are taken and then the parameters are estimated for each re-sample. The function allows for parallel calculations, via mclapply(qv), which reduces some of the computational burden. To use parallel computing, specify mc.cores>1. Note that this will not work on Windows computers, as mclapply(qv) will not work.
#' The Bayesian bootstrap method is operationally equivalent to the simple bootstrap, except that the weightings are non-integral. See Rubin (1981). It might be tempting to reduce the tolerance for convergence of each estimation procedure. We recommend not doing this as it is likely to have the effect of artificially reducing estimates of uncertainty. This occurs as the resampled estimates are likely to be closer to their starting values (the MLEs from the original data set).
#' @return An object of class "species_mix.bootstrap", which is essentially a matrix with nboot rows and the number of columns equal to the number of parameters matrix. Each row gives a bootstrap estimate of the parameters.
#' @references Rubin, D.B. (1981) The Bayesian Bootstrap. The Annals of Statistics \emph{9}:130--134.
#' @importFrom stats vcov
#' @rdname bootstrap
#' @export bootstrap
"bootstrap" <- function (object, nboot=10, type="BayesBoot", mc.cores=1, quiet=TRUE, ...){
  UseMethod("bootstrap", object)
}

#' @export
"bootstrap.species_mix" <-function (object, nboot=10, type="BayesBoot",
                                    mc.cores=1, quiet=FALSE,...){
  if (nboot < 1)
    stop( "No Bootstrap samples requested.  Please set nboot to something > 1.")
  if( ! type %in% c("BayesBoot","SimpleBoot"))
      stop( "Unknown Bootstrap type, choices are BayesBoot and SimpleBoot.")

  n.reorder <- 0
  object$titbits$control$optimise <- TRUE #just in case it was turned off

  if( !quiet){
      chars <- c("><(('> ","_@_'' ")
      pb <- txtProgressBar(min = 1, max = nboot, style = 3, char = chars[sample(length(chars),1)])
  }

   if( type == "SimpleBoot"){
      all.wts <- matrix( sample( 1:object$S, nboot*object$S, replace=TRUE), nrow=nboot, ncol=object$S)
      tmp <- apply( all.wts, 1, table)
      all.wts <- matrix( 0, nrow=nboot, ncol=object$S)
      for( ii in seq_along( tmp))
        all.wts[ii, as.numeric( names( tmp[[ii]]))] <- tmp[[ii]]
    }
    if( type == "BayesBoot")
      all.wts <- object$S * gtools::rdirichlet( nboot, rep( 1, object$S))
    my.inits <- setup_inits_sam(object$coefs, object$S, object$G,
                                object$titbits$X, object$titbits$W,
                                object$titbits$U, object$disty,
                                return_list = TRUE)

    tmpOldQuiet <- object$titbits$control$quiet
    object$titbits$control$quiet <- TRUE

    my.fun <- function(dummy){
      if( !quiet) setTxtProgressBar(pb, dummy)
      # disty_cases <- c("bernoulli","binomial","poisson", "ippm", "negative.binomial", "tweedie", "gaussian")
      disty <- object$disty
      linky <- object$link
      dumbOut <- capture.output(
        samp.object <- ecomix::species_mix.fit(y=object$titbits$Y,
                                               X=object$titbits$X,
                                               W=object$titbits$W,
                                               U=object$titbits$U,
                                               offset = object$titbits$offset,
                                               spp_weights = all.wts[dummy,,drop=TRUE],
                                               site_spp_weights = object$titbits$site_spp_weights,
                                               G = object$G,
                                               S = object$S,
                                               # y_is_na = object$titbits$y_is_na,
                                               disty = disty,
                                               linky = linky,
                                               size = object$titbits$size,
                                               powers = object$titbits$powers,
                                               control = object$titbits$control,
                                               inits = my.inits)
        )
      boot.res <- c(samp.object$alpha,samp.object$beta,samp.object$eta,
                    samp.object$gamma,samp.object$delta,samp.object$theta)
      return(boot.res)
    }

    tmp <- plapply(seq_len(nboot), my.fun, .parallel = mc.cores)
    boot.estis <- do.call( "rbind", tmp)
  # }

  boot.estis <- boot.estis[,!apply(boot.estis=='-999999', 2, any,na.rm=TRUE)]
  object$titbits$control$quiet <- tmpOldQuiet
  class( boot.estis) <- "species_mix.bootstrap"
  return( boot.estis)
}

"species_mix_boot_parametric" <- function( object, nboot){
  if( nboot > 0){
    if( is.null( object$vcov)){
      message( "An estimate of the variance matrix for regression parameters is required. Please run object$vcov <- vcov(), see ?vcov.species_mix for help")
      return( NULL)
    }
    allCoBoot <- my.rmvnorm( n=nboot, mean=as.numeric(unlist( object$coefs)), sigma=object$vcov, method='eigen')
    return( allCoBoot)
  }
  else{
    boot.estis <- matrix( unlist( object$coef), nrow=1)
    return( boot.estis)
  }
}
