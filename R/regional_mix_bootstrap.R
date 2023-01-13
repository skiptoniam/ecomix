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

