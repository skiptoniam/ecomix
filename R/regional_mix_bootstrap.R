#' @rdname bootstrap.regional_mix
#' @name bootstrap.regional_mix
#' @title Performs bootstrap sample and estimation for regimix objects. Useful for calculating measures of uncertainty in predictions from a regimix object, and also about the regimix parameter estimates. This function can be used in conjunction with vcov.regimix(qv) and predict.regimix(qv). In partciular, these bootstrap samples can be used to gauge variability in parameter estimates and hence the model itself.
#' @aliases regional_mix.bootstrap
#' @title Bootstraps a regional_mix object.
#' @description Performs bootstrap sample and estimation for regional_mix objects. Useful for calculating measures of uncertainty in predictions from a regional_mix object, and also about the regional_mix parameter estimates. This function can be used in conjunction with vcov.regional_mix(qv) and predict.regional_mix(qv). In partciular, these bootstrap samples can be used to gauge variability in parameter estimates and hence the model itself.
#' @param object an object obtained from fitting a RCP mixture model (class "regional_mix"). Such as that generated from a call to regional_mix(qv). This object will need to be created with the argument titbits=TRUE as these pieces of data are needed in the re-fitting. Of course, you could try to just save parts of titbits but the memory-hungry data is all required.
#' @param nboot a numeric scalar giving the number of bootstrap samples to obtain. More is better, but takes longer.
#' @param type a character string giving the type of bootstrap to perform. Options are:"SimpleBoot" which gives sample resampling, and "BayesBoot" (default) which gives Bayesian Bootstrap sampling. The nomenclature, and the Bayesian Bootstrap, come from Rubin (1981).
#' @param mc.cores an integer giving the number of cores to run the bootstrap samples on. The default is 1, that is no parallelisation. This parameter is redundant on Windows machines as the method of parallelisation, mclapply(), is not available there.
#' @param quiet should the progress bar be printed to the output device?If quiet=FALSE (default) then the progress bar is printed.
#' @param MLstart should each bootstrap estimation start at the original model's ML estimate? Default is TRUE for \code{yes} it should.
#' @param \\dots Ignored
#' @description This function can take a while to run -- it is a bootstrap function. nboot re-samples of the data are taken and then the parameters are estimated for each re-sample. The function allows for parallel calculations, via mclapply(qv), which reduces some of the computational burden. To use parallel computing, specify mc.cores>1. Note that this will not work on Windows computers, as mclapply(qv) will not work.
#' The Bayesian bootstrap method is operationally equivalent to the simple bootstrap, except that the weightings are non-integral. See Rubin (1981). It might be tempting to reduce the tolerance for convergence of each estimation procedure. We recommend not doing this as it is likely to have the effect of artificially reducing estimates of uncertainty. This occurs as the resampled estimates are liekly to be closer to their starting values (the MLEs from the original data set).
#' @return An object of class "regional_mix.bootstrap", which is essentially a matrix with nboot rows and the number of columns equal to the number of parameters matrix. Each row gives a bootstrap estimate of the parameters.
#' @references Foster, S.D., Lyons, M. and Hill, N. (in prep.) Ecological Groupings of Sample Sites in the presence of sampling artefacts.
#' Rubin, D.B. (1981) The Bayesian Bootstrap. The Annals of Statistics \emph{9}:130--134.
#' @export

"bootstrap.regional_mix" <-function (object,
                                     nboot=10,
                                     type="BayesBoot",
                                     mc.cores=1,
                                     quiet=FALSE,
                                     # orderSamps=FALSE,
                                     MLstart=TRUE,
                                     ...){
  if (nboot < 1)
    stop( "No Boostrap samples requested.  Please set nboot to something > 1.")
  if( ! type %in% c("BayesBoot","SimpleBoot"))
    stop( "Unknown boostrap type, choices are BayesBoot and SimpleBoot.")
  n.reorder <- 0
  object$titbits$control$optimise <- TRUE #just in case it was turned off (see regional_mix.multfit)
  if( !quiet){
    chars <- c("><(('> ","_@_'' ")
    pb <- txtProgressBar(min = 1, max = nboot, style = 3, char = chars[sample(length(chars),1)]) }
  if( type == "SimpleBoot"){
    all.wts <- matrix( sample( 1:object$n, nboot*object$n, replace=TRUE), nrow=nboot, ncol=object$n)
    tmp <- apply( all.wts, 1, table)
    all.wts <- matrix( 0, nrow=nboot, ncol=object$n)
    for( ii in seq_along( tmp))
      all.wts[ii, as.numeric( names( tmp[[ii]]))] <- tmp[[ii]]
  }
  if( type == "BayesBoot")
    all.wts <- object$n * gtools::rdirichlet( nboot, rep( 1, object$n))
  if( MLstart)
    my.inits <- unlist( object$coef)
  else{
    my.inits <- "random"
    orderSamps <- TRUE
  }

  my.fun <- function(dummy){
    if( !quiet)
      setTxtProgressBar(pb, dummy)
    dumbOut <- capture.output(
      samp.object <- regional_mix.fit( outcomes=object$titbits$Y, W=object$titbits$W, X=object$titbits$X, offy=object$titbits$offset, wts=object$titbits$wts * all.wts[dummy,,drop=TRUE], disty=object$titbits$disty, nRCP=object$nRCP, power=object$titbits$power, inits=my.inits, control=object$titbits$control, n=object$n, S=object$S, p.x=object$p.x, p.w=object$p.w))
    return( unlist( samp.object$coef))
  }

  flag <- TRUE
  tmpOldQuiet <- object$titbits$control$quiet
  object$titbits$control$quiet <- TRUE
  if( Sys.info()['sysname'] == "Windows" | mc.cores==1){
    boot.estis <- matrix(NA, nrow = nboot, ncol = length(unlist(object$coef)))
    for (ii in 1:nboot) {
      if( !quiet)
        setTxtProgressBar(pb, ii)
      boot.estis[ii, ] <- my.fun( ii)
    }
    flag <- FALSE
  }
  if( flag){  #has this already been done sequentially?
    if( !quiet)
      message( "Progress bar may not be monotonic due to the vaguaries of parallelisation")
    tmp <- parallel::mclapply( 1:nboot, my.fun, mc.silent=quiet, mc.cores=mc.cores)
    boot.estis <- do.call( "rbind", tmp)
  }
  object$titbits$control$quiet <- tmpOldQuiet
  if( !quiet)
    message( "")
  colnames( boot.estis) <- get_long_names_rcp( object)
  class( boot.estis) <- "regional_mix.bootstrap"
  return( boot.estis)
}

"regional_mix_bootParametric" <- function( fm, mf, nboot){
  if( nboot > 0){
    if( is.null( fm$vcov)){
      message( "An estimate of the variance matrix for regression parameters is required. Please run fm$vcov <- vcov(), see ?vcov.regional_mix for help")
      return( NULL)
    }
    allCoBoot <- my.rmvnorm( n=nboot, mean=as.numeric( unlist( fm$coefs)), sigma=fm$vcov, method='eigen')
    return( allCoBoot)
  }
  else{
    boot.estis <- matrix( unlist( fm$coef), nrow=1)
    return( boot.estis)
  }
}

