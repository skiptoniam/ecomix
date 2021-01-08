##### S3 Class exports #####
#' @rdname AIC.regional_mix
#' @name AIC.regional_mix
#' @title AIC from a regional_mix model
#' @param object A regional_mix object
#' @param \\dots extra things for AIC
#' @param k AIC penality
#' @export
"AIC.regional_mix" <- function (object, ..., k = 2){
  p <- length(unlist(object$coefs))
  if (is.null(k))
    k <- 2
  star.ic <- -2 * object$logl + k * p
  return(star.ic)
}

#' @rdname BIC.regional_mix
#' @name BIC.regional_mix
#' @title BIC from a regional_mix model
#' @param object A regional_mix object
#' @param \\dots Ignored
#' @export
#' @importFrom stats BIC
"BIC.regional_mix" <- function (object, ...){
  p <- length(unlist(object$coefs))
  k <- log(object$n)
  star.ic <- -2 * object$logl + k * p
  return(star.ic)
}

#' @rdname coef.regional_mix
#' @name coef.regional_mix
#' @title Print the coefficientf from a regional_mix model
#' @param object A regional_mix object
#' @param \\dots Ignored
#'@export

"coef.regional_mix" <- function (object, ...){
  res <- list()
  res$alpha <- object$coefs$alpha
  names( res$alpha) <- object$names$spp
  if( !is.null( object$coef$tau)){
    res$tau <- matrix(object$coefs$tau, nrow = object$nRCP - 1, ncol = object$S)
    colnames( res$tau) <- object$names$spp
  }
  if( !is.null( object$coef$beta)){
    res$beta <- matrix(object$coefs$beta, nrow = object$nRCP - 1, ncol = object$p.x)
    colnames( res$beta) <- object$names$Xvars
  }
  if( !any(is.null( object$coef$gamma),length( object$coef$gamma)==0)){
    res$gamma <- matrix( object$coef$gamma, nrow=object$S, ncol=object$p.w)
    colnames( res$gamma) <- object$names$Wvars
    rownames( res$gamma) <- object$names$spp
  }
  if( !is.null( object$coef$disp)){
    res$logDisp <- object$coef$disp
    names( res$logDisp) <- object$names$spp
  }

  return(res)
}

#'@title cooks.distance
#'@rdname cooks.distance
#'@description Performs leave-some-out measures for a regional_mix model. This includes a measure of how much effect leaving out an observation has on the probability of each site's RCP label. Also, this function can be used as a cross-validation workhorse.

#'@param model A regional_mix object whose fit you want to assess
#'@param oosSize The size of the with held partitions (out-of-sample size). Use 1 (default) for leave-one-out statistics, such as Cook's distance and leave-one-out validation.
#'@param times The number of times to perform the re-estimation (the number of leave out groups). For each 1:times a random partition of the data, of size oosSize, is taken and the model is fitted to one of the partitions. It is predicted to the other partition. The exception is when oosSize=1 and times=model$n (leave-one-out). In such cases (the default too), the observations are left out one-by-one and not randomly.
#'@param mc.cores The number of cores to spread the workload over. Default is 1. Argument is useless on Windows machines ??? see ?parallel::mclapply
#'@param quiet Should printing be suppressed? Default is no, it should not. Note that in either case, printing of the iteration trace etc is suppressed for each regional_mix fit.
#'@param \\dots ignored
#'@return An object of class regiCooksD. It is a list of 4 elements:
#'@return Y the species data,
#'@return CV the model$n by model$S by times array of out-of-sample predictions (this array contains a lot of NAs for where predictions would in-sample),
#'@return cooksD a model$n by model$nRCP matrix of statistics that resemble Cook's distance. The statistic is the change in the prediction of RCP probability from the model with all the data to the model with only the in-sample data, and predLogL the predictive log-likelihood of each point in each withheld sample (log-likelihood contributions of withheld observations, again there will be many NAs).
#'@export
#'@examples
#' \dontrun{
#' #not run as R CMD check complains about the time taken.
#' #This code will take a little while to run (<1 minute on my computer)
#' #For leave-one-out cooks distance, use oosSize=1
#' #for serious use times will need to be larger.
#' system.time({
#'   example( regional_mix);
#'   cooksD <- cooks.distance( fm, oosSize=10, times=25)
#' })
#' #For leave-one-out cooks distance, use oosSize=1
#' cooksD <- cooks.distance( fm, oosSize=10, times=5)
#' }


"cooks.distance.regional_mix" <- function( model, ..., oosSize=1, times=model$n, mc.cores=1, quiet=FALSE){
  if (oosSize > model$n %/% 2)
    stop("Out of sample is more than half the size of the data! This is almost certainly an error.  Please set `oosSize' to something smaller.")
  if (is.null(model$titbits))
    stop("Model doesn't contain all information required for cross validation.  Please supply model with titbits (from titbits=TRUE in regional_mix call)")
  if ( !quiet)
    pb <- txtProgressBar(min = 1, max = times, style = 3, char = "><(('> ")

  funny <- function(x) {
    if (!quiet)
      setTxtProgressBar(pb, x)
    if( oosSize!=1 | times!=model$n) #do we need to sample?
      OOBag <- sample(1:model$n, oosSize, replace = FALSE)
    else
      OOBag <- x
    inBag <- (1:model$n)[!(1:model$n) %in% OOBag]
    new.wts <- model$titbits$wts
    new.wts[OOBag] <- 0
    control <- model$titbits$control
    control$quiet <- TRUE
    control$trace <- 0
    control$optimise <- TRUE
    tmpmodel <- regional_mix.fit(outcomes = model$titbits$Y,
                                 W = model$titbits$W, X = model$titbits$X, offy = model$titbits$offset,
                                 wts = new.wts, disty = model$titbits$disty, nRCP = model$nRCP,
                                 power = model$titbits$power, inits = unlist(model$coef),
                                 control = control, n = model$n, S = model$S, p.x = model$p.x,
                                 p.w = model$p.w)
    OOSppPreds <- matrix(NA, nrow = tmpmodel$n, ncol = tmpmodel$S)
    for (ss in 1:tmpmodel$S)
      OOSppPreds[OOBag, ss] <- rowSums(tmpmodel$mus[OOBag, ss,] * tmpmodel$pis[OOBag, , drop=FALSE])
    newPis <- tmpmodel$pis
    r.negi <- model$pis - newPis
    r.negi[OOBag,] <- NA
    r.negi <- colMeans( r.negi, na.rm=TRUE)
    #great lengths to calc pred logl...
    #great lengths indeed...
    alpha.score <- as.numeric(rep(NA, model$S))
    tau.score <- as.numeric(matrix(NA, ncol = model$S, nrow = model$nRCP - 1))
    beta.score <- as.numeric(matrix(NA, ncol = ncol(model$titbits$X), nrow = model$nRCP - 1))
    if( model$p.w > 0){
      gamma.score <- as.numeric(matrix( NA, nrow=model$S, ncol=model$p.w))
      gamma <- tmpmodel$coef$gamma
      W <- model$titbits$W
    }
    else
      gamma.score <- W <- gamma <- -999999
    if( model$titbits$disty %in% 3:5){
      disp.score <- as.numeric( rep( NA, model$S))
      disp <- stats::coef(model)$logDisp
    }
    else
      disp.score <- -999999
    scoreContri <- -999999
    #model quantities
    #    pis <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))  #container for the fitted RCP model
    #    mus <- as.numeric(array( NA, dim=c( n, S, nRCP)))  #container for the fitted spp model
    logCondDens <- as.numeric(matrix(NA, nrow = model$n, ncol = model$nRCP))
    logls <- as.numeric(rep(NA, model$n))
    conv <- as.integer(0)
    tmplogl <- .Call("RCP_C", as.numeric( model$titbits$Y), as.numeric(model$titbits$X), as.numeric( model$titbits$W), as.numeric(model$titbits$offset), as.numeric(model$titbits$wts),
                     as.integer(model$S), as.integer(model$nRCP), as.integer(model$p.x), as.integer(model$p.w), as.integer(model$n), as.integer( model$titbits$disty),
                     as.numeric( tmpmodel$coef$alpha), as.numeric( tmpmodel$coef$tau), as.numeric( tmpmodel$coef$beta), as.numeric( gamma), as.numeric( tmpmodel$coef$disp), as.numeric( model$titbits$power),
                     as.numeric(model$titbits$control$penalty), as.numeric(model$titbits$control$penalty.tau), as.numeric(model$titbits$control$penalty.gamma), as.numeric(model$titbits$control$penalty.disp[1]), as.numeric(model$titbits$control$penalty.disp[2]),
                     alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                     as.numeric( tmpmodel$pis), as.numeric( tmpmodel$mus), logCondDens, logls,
                     as.integer(model$titbits$control$maxit), as.integer(model$titbits$control$trace), as.integer(model$titbits$control$nreport), as.numeric(model$titbits$control$abstol), as.numeric(model$titbits$control$reltol), as.integer(conv),
                     as.integer(FALSE), as.integer(TRUE), as.integer(FALSE), as.integer(FALSE), as.integer(FALSE), PACKAGE = "ecomix")
    ret.logl <- rep( NA, model$n)
    ret.logl[OOBag] <- logls[OOBag]

    return( list( OOSppPreds=OOSppPreds, cooksDist=r.negi, predLogL=ret.logl))
  }
  if (!quiet & mc.cores>1 & Sys.info()['sysname'] != "Windows")
    message("Progress bar may not be monotonic due to the vaguaries of parallelisation")
  tmp <- parallel::mclapply(1:times, funny, mc.cores = mc.cores)
  if (!quiet)
    message("")
  cooksD <- t( sapply( tmp, function(x) x$cooksDist))
  OOpreds <- array(NA, dim = c(model$n, model$S, times), dimnames = list(rownames(model$titbits$X), colnames(model$titbits$Y), paste("CVset", 1:times, sep = "")))
  for (bb in 1:times)
    OOpreds[, , bb] <- tmp[[bb]]$OOSppPreds
  logls <- sapply( tmp, function(x) x$predLogL)
  colnames( logls) <- rownames( cooksD) <- paste( "OOS",1:times,sep="_")
  ret <- list(Y = model$titbits$Y, CV = OOpreds, cooksD=cooksD, predLogL=logls)
  class(ret) <- "regiCooksD"

  return(ret)
}

#' @rdname extractAIC.regional_mix
#' @name extractAIC.regional_mix
#' @title Extract AIC from regional_mix model.
#' @param fit  Fitted RCP model
#' @param scale scale parameter
#' @param k AIC parameter
#' @param \\dots Ignored
#' @export

"extractAIC.regional_mix" <- function (fit, scale = 1, k = 2, ...){
  n <- fit$n
  edf <- length(unlist(stats::coef(fit)))
  if (is.null(k))
    k <- 2
  aic <- -2 * logLik(fit) + k * edf
  return(c(edf, aic))
}


#' @rdname plot.regional_mix
#' @name plot.regional_mix
#' @title Plot residuals from a regional_mix model.
#' @param x a fitted regional_mix model you wish to plot
#' @param \\dots Additional plotting calls
#' @param type What type of residuals to plot? 'RQR'; random quantile residuals or 'deviance' residuals.
#' @param nsim Number of simulations to run
#' @param alpha.conf The bounds of the confidence intervals.
#' @param quiet Run in quiet mode.
#' @param species Which species to plot as residuals.
#' @param fitted.scale What scale to plot the residuals on?
#' @details The two types of residuals are inherently different. The "RQR" residuals produce a residual for each species at each site and the "deviance" residuals produce a site residual (no species level residual). The plots also differ, the "RQR" type generates a single normal QQ-plot for all species and all sites, and a residual versus fitted plot for all species and sites (Described in Foster et al, 2013). The "deviance" type generates a pair of Tukey mean-difference plots, similar in spirit to a QQ-plot. The first is for point-wise confidence intervals and the second is for approximate global intervals. See Foster et al (2013) for details.
#' The family for the "RQR" residuals should be standard normal. For "deviance" residuals, the distribution is unknown and simulation is used to graphically assess how odd the observed residuals look compared to ones generated assuming the model is correct.
#' @references Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.
#' Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245
#' Foster, S.D., Hill, N.A. and Lyons, M., 2017. Ecological grouping of survey sites when sampling artefacts are present. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), pp.1031-1047.
#' @export

"plot.regional_mix" <- function (x, ..., type="RQR", nsim = 100,
                                 alpha.conf = c(0.9, 0.95, 0.99), quiet=FALSE,
                                 species="AllSpecies", fitted.scale="response"){
  if( ! type %in% c("RQR","deviance"))
    stop( "Unknown type of residuals. Options are 'RQR' and 'deviance'.\n")
  if( ! all( species %in% c("AllSpecies",x$names$spp)))
    stop( "Unknown species.  Options are 'AllSpecies' or any one of the species names as supplied (and stored in x$names$spp)")

  if( type=="deviance"){
    obs.resid <- residuals( x, type="deviance")
    shad <- rev(seq(from = 0.8, to = 0.5, length = length(alpha.conf)))
    allResids <- matrix(NA, nrow = x$n, ncol = nsim)
    X <- x$titbits$X
    p.x <- ncol( X)
    if( class( x$titbits$species_formula)=="formula"){
      form.W <- x$titbits$species_formula
      W <- x$titbits$W
      p.w <- ncol( W)
    }
    else{
      form.W <- NULL
      W <- -999999
      p.w <- 0
    }
    offy <- x$titbits$offset
    wts <- x$titbits$wts
    Y <- x$titbits$Y
    disty <- x$titbits$disty
    power <- x$titbits$power
    S <- x$S
    nRCP <- x$nRCP
    p.x <- x$p.x
    p.w <- x$p.w
    n <- x$n
    disty <- x$titbits$disty
    control <- x$titbits$control
    pis <- as.numeric( matrix( -999999, nrow = n, ncol = nRCP))
    mus <- as.numeric( array( -999999, dim=c( n, S, nRCP)))
    logCondDens <- as.numeric( matrix( -999999, nrow = n, ncol = nRCP))
    logls <- as.numeric(rep(-999999, n))
    alpha.score <- as.numeric(rep(-999999, S))
    tau.score <- as.numeric(matrix(-999999, nrow = nRCP - 1, ncol = S))
    beta.score <- as.numeric(matrix(-999999, nrow = nRCP - 1, ncol = p.x))
    if( p.w > 0)
      gamma.score <- as.numeric( matrix( -999999, nrow = S, ncol = p.w))
    else
      gamma.score <- -999999
    if( !is.null( x$coef$disp))
      disp.score <- as.numeric( rep( -999999, S))
    else
      disp.score <- -999999
    conv <- FALSE
    alpha <- x$coefs$alpha
    tau <- x$coefs$tau
    beta <- x$coefs$beta
    if( !is.null( form.W))
      gamma <- x$coefs$gamma
    else
      gamma <- -999999
    if( any( !is.null( x$coef$disp)))
      disp <- x$coef$disp
    else
      disp <- -999999
    scoreContri <- as.numeric(matrix(NA, ncol = length(unlist(x$coef)), nrow = x$n))
    if( !quiet)
      pb <- txtProgressBar(min = 1, max = nsim, style = 3, char = "><(('> ")
    for (s in 1:nsim) {
      if( !quiet)
        setTxtProgressBar(pb, s)
      newy <- as.matrix( regional_mix.simulate( nRCP=nRCP, S=S, n=n, p.x=p.x, p.w=p.w, alpha=alpha, tau=tau, beta=beta, gamma=gamma, logDisps=disp, powers=power, X=X, W=W, offset=offy, family=x$family))
      tmp <- .Call("RCP_C", as.numeric(newy[, 1:S]), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
                   as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
                   alpha, tau, beta, gamma, disp, power,
                   as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
                   alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                   pis, mus, logCondDens, logls,
                   as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
                   as.integer( FALSE), as.integer( TRUE), as.integer( FALSE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")

      allResids[, s] <- get_residuals_rcp(logls, Y, x$family, x$coef, nRCP, type="deviance", powers=power, quiet=TRUE)
    }
    if( !quiet)
      message("")
    allResidsSort <- apply(allResids, 2, sort)
    quants <- c(0.5, (1 - alpha.conf)/2, alpha.conf + (1 - alpha.conf)/2)
    envel <- t(apply(allResidsSort, 1, quantile, probs = quants, na.rm = TRUE))
    sort.resid <- sort(obs.resid)
    empQuant <- envel[, 1]
    diff <- sweep(envel[, -1], 1, empQuant, "-")
    realMeans <- (sort.resid + empQuant)/2
    realDiff <- sort.resid - empQuant
    # par(mfrow = c(1, 2))
    plot(rep(realMeans, 1 + 2 * length(alpha.conf)), c(diff,
                                                       realDiff), sub = "Pointwise Confidence",
         ylab = "Observed - Expected", xlab = "(Observed+Expected)/2",
         type = "n")
    for (aa in rev(seq_along(alpha.conf))) polygon(c(realMeans, rev(realMeans)),
                                                   c(diff[, aa], rev(diff[, aa + length(alpha.conf)])),
                                                   col = grey(shad[aa]), border = NA)
    points(realMeans, realDiff, pch = 20)
    abline(h = 0)
    globEnvel <- envel
    for (ii in 2:(length(alpha.conf) + 1)) globEnvel[, ii] <- globCIFinder(x = allResidsSort, en = envel[, ii], alpha = quants[ii], nsim = nsim)
    for (ii in 1 + (length(alpha.conf) + 1):(2 * length(alpha.conf))) globEnvel[, ii] <- globCIFinder(x = allResidsSort, en = envel[, ii], alpha = quants[ii], nsim = nsim)
    empQuant <- globEnvel[, 1]
    diff <- sweep(globEnvel[, -1], 1, empQuant, "-")
    realMeans <- (sort.resid + empQuant)/2
    realDiff <- sort.resid - empQuant
    plot(rep(realMeans, 1 + 2 * length(alpha.conf)), c(diff,
                                                       realDiff), sub = "Global Confidence",
         ylab = "Observed - Expected", xlab = "(Observed+Expected)/2",
         type = "n")
    for (aa in rev(seq_along(alpha.conf)))
      polygon(c(realMeans, rev(realMeans)), c(diff[, aa], rev(diff[, aa + length(alpha.conf)])), col = grey(shad[aa]), border = NA)
    points(realMeans, realDiff, pch = 20)
    abline(h = 0)
    return(NULL)
  }

  if( type=="RQR"){
    obs.resid <- residuals( x, type="RQR", quiet=quiet)
    S <- x$S
    sppID <- rep( TRUE, S)
    if( species != "AllSpecies"){
      sppID <- x$names$spp %in% species
      obs.resid <- obs.resid[,sppID, drop=FALSE]
      S <- ncol( obs.resid)
    }
    if( sum( obs.resid==Inf | obs.resid==-Inf) > 0){
      message( "Infinite residuals removed from residual plots:", sum( obs.resid==Inf | obs.resid==-Inf), "in total.")
      obs.resid[obs.resid==Inf | obs.resid==-Inf] <- NA
    }
    spp.cols <- rep( 1:S, each=x$n)
    main <- match.call( expand.dots=TRUE)$main
    if( is.null( main)){
      if( species=="AllSpecies")
        main <- "All Residuals"
      else
        if( length( species)==1)
          main<-species
        else
          main<-""
    }
    sub <- match.call( expand.dots=TRUE)$sub
    if( is.null( sub))
      sub <- "Colours separate species"
    par( mfrow=c(1,2))
    qqnorm(obs.resid, col=spp.cols, pch=20, main=main, sub=sub)
    #    qqline( obs.resid)  #this doesn't actually poduce a y=x line.  It is only(?) appropriate if the scales of the two sets are different.
    abline( 0,1,lwd=2)

    preds <- matrix( NA, nrow=x$n, ncol=S)
    for( ii in 1:x$n){
      preds[ii,] <- rowSums( x$mu[ii,sppID,] * matrix( rep( x$pi[ii,], each=S), nrow=S, ncol=x$nRCP))
    }
    switch( fitted.scale,
            log = { loggy <- "x"},
            logit = { loggy <- ""; preds <- log( preds / (1-preds))},
            {loggy <- ""})
    plot( preds, obs.resid, xlab="Fitted", ylab="RQR", main="Residual versus Fitted", sub="Colours separate species", pch=20, col=rep( 1:S, each=x$n), log=loggy)
    abline( h=0)

  }
}

#' @rdname plot.regional_mix.sampling
#' @name plot.regional_mix.sampling
#' @title Plot residuals from a regional_mix model.
#' @param object a fitted regional_mix model you wish to plot
#' @param object2 a bootstrap object from regional_mix.bootstrap
#' @param type What type of residuals to plot? 'RQR'; random quantile residuals or 'deviance' residuals.
#' @param nsim Number of simulations to run
#' @param alpha.conf The bounds of the confidence intervals.
#' @param quiet Run in quiet mode.
#' @param species Which species to plot as residuals.
#' @param fitted.scale What scale to plot the residuals on?
#' @details The two types of residuals are inherently different. The "RQR" residuals produce a residual for each species at each site and the "deviance" residuals produce a site residual (no species level residual). The plots also differ, the "RQR" type generates a single normal QQ-plot for all species and all sites, and a residual versus fitted plot for all species and sites (Described in Foster et al, 2013). The "deviance" type generates a pair of Tukey mean-difference plots, similar in spirit to a QQ-plot. The first is for point-wise confidence intervals and the second is for approximate global intervals. See Foster et al (2013) for details.
#' The family for the "RQR" residuals should be standard normal. For "deviance" residuals, the distribution is unknown and simulation is used to graphically assess how odd the observed residuals look compared to ones generated assuming the model is correct.
#' @references Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.
#' Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245
#' Foster, S.D., Hill, N.A. and Lyons, M., 2017. Ecological grouping of survey sites when sampling artefacts are present. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), pp.1031-1047.
#' @export

"plot.regional_mix.sampling" <-function(x, ... , object2=NULL, sampling_levels,
                                        CI=c(0.025, 0.975), col="black", lty=1){

  require(lattice)
  gammas<-grepl("gamma",dimnames(object2)[[2]])
  temp_dat<-object2[,gammas]

  temp<-data.frame(avs=as.numeric(unname(colMeans(temp_dat))),
                   t(apply(temp_dat, 2, quantile, probs=CI)),
                   #sampling_var=rep(sampling_names, each=length(length(best_mod$names$spp))),
                   sampling_var=sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 3),
                   #Species=factor(rep(best_mod$names$spp,length(sampling_names))))
                   Species=factor(sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 1)))

  names(temp)[2:3]<-c("lower", "upper")
  temp$Species<-gsub("."," ", temp$Species, fixed=TRUE) #get rid of '.' in species names
  temp$Species<-as.factor(temp$Species) #convert back to factor
  temp$Species <- factor(temp$Species, levels=rev(levels(temp$Species)))

  trellis.par.set(superpose.symbol=list(pch=16,col=col, cex=1.2),
                  superpose.line=list(col="transparent"))
  dotplot(Species ~ avs, groups=sampling_var, data=temp, cols=col, lty=lty, low=temp$lower, high=temp$upper, subscript=TRUE,
          auto.key=list(space="top", columns=2, cex=1.4, text=legend_fact),
          ylab=list(cex=1.4), xlab=list("Coefficient",cex=1.4),
          scales = list(tck = c(1, 0), x=list(cex=1.2), y=list(cex=1.2)),
          prepanel = function(x, y, ...) { list(xlim=range(temp$lower, temp$upper)) },
          panel=panel.superpose,
          panel.groups=function(x, y, subscripts, group.number, cols, low, high, ...)
          {
            if(group.number==1) jiggle <- 0.1 else jiggle <- -0.1
            panel.abline(v=0, lty=2)
            panel.abline(h=1:length(best_mod$names$spp), col.line="light grey", lty=1)
            panel.dotplot(x, y+jiggle, group.number, ...)
            panel.arrows(low[subscripts], y+jiggle, high[subscripts], y+jiggle, code=3, angle=90,
                         length=0.05, col=cols[group.number], lty=lty[group.number])
            #panel.segments(temp$lower, y+jiggle, + temp$upper, y+jiggle, lty = lty, col =col, lwd=2, cex=1.2)
          })

}



#' @rdname plot.regional_mix_stab
#' @name plot.regional_mix_stab
#' @title Diagnostic plotting to see if RCP groups are stable
#' @description For increasing size of hold-out samples, cooks distance and predictive log-likelihood are plotted.
#' @param x x-axis
#' @param y y-axis
#' @param minWidth min width of cuts/binning
#' @param ncuts number of cuts/bins to make
#' @param ylimmo limit of y-axis
#' @param \\dots additional plotting calls
#' @export
#' @examples
#' \dontrun{
#'  #not run as R CMD check complains about the time taken.
#'  #This code will take a little while to run (about 3.5minutes on my computer)
#'  system.time(\{
#'  my.registab <- stability.regional_mix( fm, oosSizeRange=seq( from=1,to=fm$n\%/\%5,length=5),
#'                                      times=fm$n, mc.cores=2, doPlot=FALSE);
#'  plot( my.registab, minWidth=1, ncuts=15);
#'  \})
#'}

"plot.regional_mix_stab" <-function(x, y, minWidth=1, ncuts=111, ylimmo=NULL, ...){
  # par(mfrow = c(1, 2))
  matplot(c(0, x$oosSizeRange), rbind(0, x$disty), type = "b",
          ylab = "Distance from Full Model Predictions", xlab = "Number of Obs Removed",
          main = "Stability of Group Predictions", col = 1:x$nRCP,
          pch = as.character(1:x$nRCP), lty = 1)
  legend("topleft", bty = "n", lty = 1, pch = as.character(1:x$nRCP),
         col = 1:x$nRCP, legend = paste("RCP ", 1:x$nRCP, sep = ""))
  oosDiffs <- diff( c(0,x$oosSizeRange))
  oosWidth <- max( minWidth, min( oosDiffs)) / 2
  histy <- list()
  for( ii in seq_along( x$oosSizeRange)){
    tmp <- stats::na.exclude( as.numeric( x$predlogls[ii,,]))
    histy[[ii]] <- hist( tmp, breaks=ncuts, plot=FALSE)
  }
  max.dens <- max( sapply( sapply( histy, function(x) x$density), max))
  if( is.null( ylimmo))
    ylimmo <- range( sapply( histy, function(x) x$breaks))
  plot( 0, 0, ylab = "Pred LogL (OOS)", xlab = "Number of Obs Removed", main = "Stability of Pred Logl", xlim = c(0-oosWidth, max(x$oosSizeRange)+oosWidth), ylim=ylimmo, type = "n")
  for( ii in seq_along( x$oosSizeRange))
    for( jj in seq_along( histy[[ii]]$density))
      rect( xleft=x$oosSizeRange[ii]-oosWidth, xright=x$oosSizeRange[ii]+oosWidth, ybottom=histy[[ii]]$breaks[jj], ytop=histy[[ii]]$breaks[jj+1], col=rgb( colorRamp( c("#E6FFFF","blue"))(histy[[ii]]$density[jj]/max.dens), maxColorValue=255), border=NA)
  tmp <- stats::na.exclude( as.numeric( x$logl.sites))
  histy <- hist( tmp, breaks=ncuts, plot=FALSE)
  for( jj in seq_along( histy$density))
    rect( xleft=0-oosWidth, xright=0+oosWidth, ybottom=histy$breaks[jj], ytop=histy$breaks[jj+1], col=rgb( colorRamp( c("#FFE6FF","red"))(histy$density[jj]/max( histy$density)), maxColorValue=255), border=NA)
  lines(c(0, x$oosSizeRange), c(mean(x$logl.sites), apply(x$predlogls,
                                                          1, mean, na.rm = TRUE)), lwd = 2, col = "black")
  invisible(TRUE)
}

#' @rdname predict.regional_mix
#' @name predict.regional_mix
#' @title Predicts RCP probabilities at a series of sites. Confidence intervals are available too.
#' @param object an object obtained from fitting a RCP mixture model. Such as that generated from a call to regional_mix(qv).
#' @param object2 a regional_mix object obtained from bootstrapping the regional_mix object. Such as that generated from a call to regional_mix_boot(qv). If not supplied, then predict.regional_mix will do parametric bootstrapping (otherwise non-parametric bootstrap).
#' @param newdata a data.frame (or something that can be coerced) containing the values of the covariates where predictions are to be made. If NULL (the default) then predictions are made at the locations of the original data.
#' @param nboot the number of parametric bootstrap samples to take for the bootstrap predictions, standard errors and confidence intervals. The default is 0, that is no bootstrapping is to be done and point predictions only are given. If object2 is not NULL, then the number of bootstrap samples is taken from that object (this argument is then ignored).
#' @param alpha a numeric within [0,1] (well [0.5,1] really) indicating the specified confidence for the confidence interval. Argument is redundant if nboot == 0.
#' @param mc.cores the number of cores to spread the computations over. Ignored if running on a Windows machine.
#' @param \\dots additional predict calls	ignorned
#' @details This function implements two separate, and quite different, bootstrapping routines. The first, attributable to Foster et al (2013), which implements a parametric bootstrap, whereby parameters are drawn from their sampling distribution (defined by the ML estimates and their asymptotic vcov matrix). Yes, the vcov function needs to be run first and stored in the the regional_mix object as $vcov. Typically, the vcov matrix is obtained using numerical derivatives, which can be slow to calculate and somewhat unstable/erratic. This was the original suggestion and has been superceeded by the non-parametric bootstrap routine. This is described in Foster et al (in prep) and bootstraps the sampling site data repeatedly, and for each bootstrap sample the model is re-estimated. Variation in the bootstrap samples is carried forward to the prediction step to guage the uncertainty.
#' The parametric bootsrap implementation of this function can take a while to run ??? it is a bootstrap function. nboot samples of the parameters are taken and then used to predict at each set of covariates defined in newdata. Quantiles of the resulting sets of bootstrap predictions are then taken. It is the last step that really takes a while. The non-parametric version of this function should not take as long as the grunt work of bootstrapping is carried out in the regional_mix_boot(qv) function.
#' Note that this function is not implemented. It could be, using the parallel package, but it is currently not. The bulk of the bootstrap calculations are done in C++, which reduces the waiting time but parallelising it would be even better.
#' @return If nboot==0 then a n x H matrix of prior predictions (n=nrow(newdata), H=number of RCPs). Each row should sum to one.
#' \item{ if nboot!=0 then a list is returned. It has elements:}{}
#' \item{ ptPreds}{the n x H matrix of point predictions}
#' \item{ bootPreds}{the n x H matrix of bootstrap point predictions (mean of bootstrap samples)}
#' \item{ bootSEs}{the n x H matrix of bootstrap standard errors for predictions}
#' \item{ bootCIs}{the n x H x 2 array of bootstrap confidence intervals. Note that bootCIs[,,1] gives the lower CIs and bootCIs[,,2] gives the upper CIs.}
#' @export

"predict.regional_mix" <- function (object, object2 = NULL, ..., newdata = NULL, nboot = 0, alpha = 0.95, mc.cores = 1){
  if (is.null(newdata)) {
    X <- object$titbits$X
    if (class(object$titbits$species_formula) == "formula") {
      form.W <- object$titbits$species_formula
      W <- object$titbits$W
      p.w <- ncol(W)
    }
    else {
      form.W <- NULL
      W <- -999999
      p.w <- 0
    }
  }
  else {
    form.X <- as.formula(object$titbit$rcp_formula)
    if (length(form.X) == 3)
      form.X[[2]] <- NULL
    X <- model.matrix(form.X, stats::model.frame(form.X, data = as.data.frame(newdata)))
    if (class(object$titbits$species_formula) == "formula") {
      W <- model.matrix(object$titbits$species_formula, stats::model.frame(object$titbits$species_formula,
                                                                           data = as.data.frame(newdata)))
      p.w <- ncol(W)
    }
    else {
      form.W <- NULL
      W <- -999999
      p.w <- 0
    }
  }
  offy <- rep(0, nrow(X))
  S <- object$S
  G <- object$nRCP
  n <- nrow(X)
  p.x <- object$p.x
  p.w <- object$p.w
  if (is.null(object2)) {
    if (nboot > 0) {
      if( !object$titbits$control$quiet)
        message("Using a parametric bootstrap based on the ML estimates and their vcov")
      my.nboot <- nboot
    }
    else
      my.nboot <- 0
    allCoBoot <- regional_mix_bootParametric(fm = object, mf = mf,
                                             nboot = my.nboot)
  }
  else {
    if( !object$titbits$control$quiet)
      message("Using supplied regional_mix_boot object (non-parametric bootstrap)")
    allCoBoot <- as.matrix(object2)
    nboot <- nrow(object2)
  }
  if (is.null(allCoBoot))
    return(NULL)
  alphaBoot <- allCoBoot[, 1:S,drop=FALSE]
  tauBoot <- allCoBoot[, S + 1:((G - 1) * S),drop=FALSE]
  betaBoot <- allCoBoot[, S + (G - 1) * S + 1:((G - 1) * p.x),drop=FALSE]
  alphaIn <- c(NA, as.numeric(object$coefs$alpha))
  alphaIn <- alphaIn[-1]
  tauIn <- c(NA, as.numeric(object$coef$tau))
  tauIn <- tauIn[-1]
  betaIn <- c(NA, as.numeric(object$coef$beta))
  betaIn <- betaIn[-1]
  if (class(object$titbits$species_formula) == "formula") {
    gammaIn <- c(NA, as.numeric(object$coef$gamma))
    gammaIn <- gammaIn[-1]
  }
  else gammaIn <- -999999
  if (any(!is.null(object$coef$disp))) {
    dispIn <- c(NA, as.numeric(object$coef$disp))
    dispIn <- dispIn[-1]
  }
  else dispIn <- -999999
  powerIn <- c(NA, as.numeric(object$titbits$power))
  powerIn <- powerIn[-1]
  predCol <- G
  ptPreds <- as.numeric(matrix(NA, nrow = n, ncol = predCol))
  bootPreds <- as.numeric(array(NA, c(n, predCol, nboot)))
  conc <- as.numeric(NA)
  mysd <- as.numeric(NA)
  outcomes <- matrix(NA, nrow = nrow(X), ncol = S)
  myContr <- object$titbits$control
  nam <- paste("RCP", 1:G, sep = "_")
  boot.funny <- function(seg) {
    if (any(segments <= 0)) {
      nboot <- 0
      bootSampsToUse <- 1
    }
    else {
      nboot <- segments[seg]
      bootSampsToUse <- (sum( segments[1:seg])-segments[seg]+1):sum(segments[1:seg])
    }
    bootPreds <- as.numeric(array(NA, c(n, predCol, nboot)))
    tmp <- .Call("RCP_predict_C", as.numeric(-999999), as.numeric(X),
                 as.numeric(W), as.numeric(offy), as.numeric(object$titbits$wts),
                 as.integer(S), as.integer(G), as.integer(p.x), as.integer(p.w),
                 as.integer(n), as.integer(object$titbits$disty),
                 as.numeric(alphaIn), as.numeric(tauIn), as.numeric(betaIn),
                 as.numeric(gammaIn), as.numeric(dispIn), as.numeric(powerIn),
                 as.numeric(myContr$penalty), as.numeric(myContr$penalty.tau),
                 as.numeric(myContr$penalty.gamma), as.numeric(myContr$penalty.disp[1]),
                 as.numeric(myContr$penalty.disp[2]), as.numeric(alphaBoot[bootSampsToUse,]),
                 as.numeric(tauBoot[bootSampsToUse,]), as.numeric(betaBoot[bootSampsToUse,]),
                 as.integer(nboot), as.numeric(ptPreds), as.numeric(bootPreds), as.integer(1),
                 PACKAGE = "ecomix")
    if (nboot == 0) {
      ret <- matrix(ptPreds, nrow = nrow(X), ncol = predCol)
      colnames(ret) <- nam
      return(ret)
    }
    bootPreds <- matrix(bootPreds, nrow = nrow(X) * predCol,
                        ncol = nboot)
    return(bootPreds)
  }
  segments <- -999999
  ret <- list()
  ptPreds <- boot.funny(1)
  if (nboot > 0) {
    if (Sys.info()["sysname"] == "Windows") {
      if( !object$titbits$control$quiet)
        message("Parallelised version of function not available for Windows machines. Reverting to single processor.")
      mc.cores <- 1
    }
    segments <- rep(nboot%/%mc.cores, mc.cores)
    if( nboot %% mc.cores > 0)
      segments[1:(nboot%%mc.cores)] <- segments[1:(nboot%%mc.cores)] + 1

    tmp <- parallel::mclapply(1:mc.cores, boot.funny, mc.cores = mc.cores)
    bootPreds <- do.call("cbind", tmp)
    bPreds <- list()
    row.exp <- rowMeans(bootPreds)
    tmp <- matrix(row.exp, nrow = nrow(X), ncol = predCol)
    bPreds$fit <- tmp
    tmp <- sweep(bootPreds, 1, row.exp, "-")
    tmp <- tmp^2
    tmp <- sqrt(rowSums(tmp)/(nboot - 1))
    tmp <- matrix(tmp, nrow = nrow(X), ncol = predCol)
    bPreds$ses <- tmp
    colnames(bPreds$fit) <- colnames(bPreds$ses) <- nam
    tmp.fun <- function(x) return(quantile(bootPreds[x, ],
                                           probs = c(0, alpha) + (1 - alpha)/2, na.rm = TRUE))
    tmp1 <- parallel::mclapply(seq_len(nrow(bootPreds)), tmp.fun,
                               mc.cores = mc.cores)
    tmp1 <- do.call("rbind", tmp1)
    tmp1 <- array(tmp1, c(nrow(X), predCol, 2), dimnames = list(NULL,
                                                                NULL, NULL))
    bPreds$cis <- tmp1[, 1:predCol, ]
    dimnames(bPreds$cis) <- list(NULL, nam, c("lower", "upper"))
    ret <- list(ptPreds = ptPreds, bootPreds = bPreds$fit,
                bootSEs = bPreds$ses, bootCIs = bPreds$cis)
  }
  else ret <- ptPreds
  gc()
  return(ret)
}

#'@rdname print.regional_mix
#'@name print.regional_mix
#'@title Prints some attributes of a regimix object.
#'@param x A fitted regional_mix object
#'@param \\dots Ignored
#'@description A list is returned that will be printed on exit, if not assigned to anything. It contains the function call and the estimated coefficients.
#'@export
#'
"print.regional_mix" <- function (x, ...){
  ret <- list()
  ret$Call <- x$call
  ret$Distribution <- x$family
  ret$coef <- stats::coef(x)
  print( ret)
  invisible(ret)
}

#' @rdname regional_mix_boot
#' @name regional_mix_boot
#' @title Performs bootstrap sample and estimation for regimix objects. Useful for calculating measures of uncertainty in predictions from a regimix object, and also about the regimix parameter estimates. This function can be used in conjunction with vcov.regimix(qv) and predict.regimix(qv). In partciular, these bootstrap samples can be used to gauge variability in parameter estimates and hence the model itself.
#' @aliases regional_mix_boot
#' @title Bootstraps a regional_mix object.
#' @description Performs bootstrap sample and estimation for regional_mix objects. Useful for calculating measures of uncertainty in predictions from a regional_mix object, and also about the regional_mix parameter estimates. This function can be used in conjunction with vcov.regional_mix(qv) and predict.regional_mix(qv). In partciular, these bootstrap samples can be used to gauge variability in parameter estimates and hence the model itself.
#' @param object an object obtained from fitting a RCP mixture model (class "regional_mix"). Such as that generated from a call to regional_mix(qv). This object will need to be created with the argument titbits=TRUE as these pieces of data are needed in the re-fitting. Of course, you could try to just save parts of titbits but the memory-hungry data is all required.
#' @param nboot a numeric scalar giving the number of bootstrap samples to obtain. More is better, but takes longer.
#' @param type a character string giving the type of bootstrap to perform. Options are:"SimpleBoot" which gives sample resampling, and "BayesBoot" (default) which gives Bayesian Bootstrap sampling. The nomenclature, and the Bayesian Bootstrap, come from Rubin (1981).
#' @param mc.cores an integer giving the number of cores to run the bootstrap samples on. The default is 1, that is no parallelisation. This parameter is redundant on Windows machines as the method of parallelisation, mclapply(), is not available there.
#' @param quiet should the progress bar be printed to the output device?If quiet=FALSE (default) then the progress bar is printed.
#' @param MLstart should each bootstrap estimation start at the original model's ML estimate? Default is TRUE for \code{yes} it should.
#' @description This function can take a while to run -- it is a bootstrap function. nboot re-samples of the data are taken and then the parameters are estimated for each re-sample. The function allows for parallel calculations, via mclapply(qv), which reduces some of the computational burden. To use parallel computing, specify mc.cores>1. Note that this will not work on Windows computers, as mclapply(qv) will not work.
#' The Bayesian bootstrap method is operationally equivalent to the simple bootstrap, except that the weightings are non-integral. See Rubin (1981). It might be tempting to reduce the tolerance for convergence of each estimation procedure. We recommend not doing this as it is likely to have the effect of artificially reducing estimates of uncertainty. This occurs as the resampled estimates are liekly to be closer to their starting values (the MLEs from the original data set).
#' @return An object of class "regional_mix_boot", which is essentially a matrix with nboot rows and the number of columns equal to the number of parameters matrix. Each row gives a bootstrap estimate of the parameters.
#' @references Foster, S.D., Lyons, M. and Hill, N. (in prep.) Ecological Groupings of Sample Sites in the presence of sampling artefacts.
#' Rubin, D.B. (1981) The Bayesian Bootstrap. The Annals of Statistics \emph{9}:130--134.
#' @export

"regional_mix_boot" <-function (object,
                                nboot=1000,
                                type="BayesBoot",
                                mc.cores=1,
                                quiet=FALSE,
                                # orderSamps=FALSE,
                                MLstart=TRUE){
  if (nboot < 1)
    stop( "No Boostrap samples requested.  Please set nboot to something > 1.")
  if( ! type %in% c("BayesBoot","SimpleBoot"))
    stop( "Unknown boostrap type, choices are BayesBoot and SimpleBoot.")
  n.reorder <- 0
  object$titbits$control$optimise <- TRUE #just in case it was turned off (see regional_mix.multfit)
  #  object$titbits$control$reltol <- max(1e-05, object$titbits$control$reltol)
  #  if( object$p.w>0)
  #    orig.data <- data.frame( cbind( object$titbits$Y, object$titbits$X, object$titbits$W, offset=object$titbits$offset, weights=rep(0,nrow(object$titbits$Y))))
  #  else
  #    orig.data <- data.frame( cbind( object$titbits$Y, object$titbits$X, offset=object$titbits$offset), weights=rep(0,nrow(object$titbits$Y)))
  if( !quiet){
    chars <- c("><(('> ","_@_'' ","@(*O*)@ ")
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
    # if( orderSamps)
    #   samp.object <- orderPost( samp.object, object)
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
  if( flag){  #has this already been done sequencially?
    if( !quiet)
      message( "Progress bar may not be monotonic due to the vaguaries of parallelisation")
    tmp <- parallel::mclapply( 1:nboot, my.fun, mc.silent=quiet, mc.cores=mc.cores)
    #    if( !quiet)
    #      message("")
    boot.estis <- do.call( "rbind", tmp)
  }
  object$titbits$control$quiet <- tmpOldQuiet
  if( !quiet)
    message( "")
  colnames( boot.estis) <- get_long_names_rcp( object)
  class( boot.estis) <- "regional_mix_boot"
  return( boot.estis)
}

#' @rdname residuals.regional_mix
#' @name residuals.regional_mix
#' @title Residuals for a regional_mix object
#' @param object A regional_mix model
#' @param \\dots Additional arguments for residuals function
#' @param type What type of residuals to plot? 'RQR'; random quantile residuals or 'deviance' residuals.
#' @param quiet Run in quiet mode.
#'@param mc.cores the number of cores to farm the jobs out to.
#' @description  The randomised quantile residuals ("RQR", from Dunn and Smyth, 1996) are defined by their marginal distribution function (marginality is over
#' other species observations within that site; see Foster et al, in prep). The result is one residual per species per site and they all should be standard
#' normal variates. Within a site they are likely to be correlated (as they share a common latent factor), but across sampling locations they will be independent.
#' The deviance residuals (as used here), are actually just square root of minus two times the log-likelihood contribution for each sampling location. We do
#' not subtract the log-likelihood of the saturated model as, at the time of writing, we are unsure what this log-likelihood should be (latent factors confuse
#' things here). This implies that the residuals will not have mean zero and their variance might also be heteroskedastic. This was not realised when writing
#' the original RCP paper (Foster et al, 2013), obviously. We still believe that these residuals have some utility, but we are unsure where that utility stops.
#' For general useage, the "RQR" residuals should probably be preferred.
#' @return
#' \item{ For type=="RQR", a number-of-sites by number-of-species matrix with the randomised quantile residuals, which should be distributed as a standard normal variate.}{}
#' \item{ For type=="deviance" a numeric vector of size object$n containing the deviance residuals.}{}
#' @references Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.
#' Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R{}. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245
#' Foster, S.D., Hill, N.A. and Lyons, M., 2017. Ecological grouping of survey sites when sampling artefacts are present. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), pp.1031-1047.
#' @export

"residuals.regional_mix" <- function( object, ..., type="RQR", quiet=FALSE, mc.cores=1){
  if( ! type %in% c("deviance","RQR","RQR.sim"))
    stop( "Unknown type of residual requested. Only deviance and RQR (for randomised quantile residuals) are implemented\n")
  if( type=="deviance"){
    resids <- sqrt( -2*object$logl.sites)
    if( !quiet){
      message( "The sign of the deviance residuals is unknown -- what does sign mean for multiple species? Their mean is also unknown -- what is a saturated model in a mixture model?")
      message( "This is not a problem if you are just looking for an over-all fit diagnostic using simulation envelopes (cf normal and half normal plots).")
      message( "It is a problem however, when you try to see how residuals vary with covariates etc.. but the meaning of these plots needs to be considered carefully as the residuals are for multiple species anyway.")
    }
  }
  if( type=="RQR"){
    resids <- matrix( NA, nrow=object$n, ncol=object$S)
    switch( object$family,
            bernoulli = { fn <- function(y,mu,logdisp,power) pbinom( q=y, size=1, prob=mu, lower.tail=TRUE)},
            poisson = { fn <- function(y,mu,logdisp,power) ppois( q=y, lambda=mu, lower.tail=TRUE)},
            negative.binomial = { fn <- function(y,mu,logdisp,power) pnbinom( q=y, mu=mu, size=1/exp( logdisp), lower.tail=TRUE)},
            tweedie = { fn <- function(y,mu,logdisp,power) fishMod::pTweedie( q=y, mu=mu, phi=exp( logdisp), p=power)},#CHECK!!!
            gaussian = { fn <- function(y,mu,logdisp,power) pnorm( q=y, mean=mu, sd=exp( logdisp), lower.tail=TRUE)})

    for( ss in 1:object$S){
      if( all( object$titbits$power==-999999))  tmpPow <- NULL else tmpPow <- object$titbits$power[ss]
      if( object$family %in% c("bernoulli","poisson","negative.binomial")){
        tmpLower <- fn( object$titbits$Y[,ss]-1, object$mus[,ss,], object$coef$disp[ss], tmpPow)
        tmpUpper <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$disp[ss], tmpPow)
        tmpLower <- rowSums( tmpLower * object$pis)
        tmpLower <- ifelse( tmpLower<0, 0, tmpLower) #get rid of numerical errors for really small negative values
        tmpLower <- ifelse( tmpLower>1, 1, tmpLower) #get rid of numerical errors for 1+epsilon.
        tmpUpper <- rowSums( tmpUpper * object$pis)
        tmpUpper <- ifelse( tmpUpper<0, 0, tmpUpper) #get rid of numerical errors for really small negative values
        tmpUpper <- ifelse( tmpUpper>1, 1, tmpUpper) #get rid of numerical errors for 1+epsilon.
        resids[,ss] <- runif( object$n, min=tmpLower, max=tmpUpper)
        resids[,ss] <- qnorm( resids[,ss])
      }
      if( object$family == "tweedie"){
        nonzero <- object$titbits$Y[,ss]>0
        tmpObs <- matrix( rep( object$titbits$Y[,ss], object$nRCP), ncol=object$nRCP)
        tmp <- matrix( fn( as.numeric( tmpObs[nonzero,]), as.numeric( object$mus[nonzero,ss,]), object$coefs$disp[ss], object$titbits$power[ss]), ncol=object$nRCP)
        tmp <- rowSums( tmp * object$pis[nonzero,])
        resids[nonzero,ss] <- qnorm( tmp)
        tmp <- matrix( fn( as.numeric( tmpObs[!nonzero,]), as.numeric( object$mus[!nonzero,ss,]), object$coefs$disp[ss], object$titbits$power[ss]), ncol=object$nRCP)
        tmp <- rowSums( tmp * object$pis[!nonzero,])
        resids[!nonzero,ss] <- qnorm( runif( sum( !nonzero), min=0, max=tmp))
      }
      if( object$family == "gaussian"){
        tmp <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$disp[ss], object$titbits$power[ss])
        tmp <- rowSums( tmp * object$pis)
        resids[,ss] <- qnorm( tmp)
      }
    }
    if( !quiet & sum( resids==Inf | resids==-Inf)>0)
      message( "Some residuals, well",sum( resids==Inf | resids==-Inf), "to be precise, are very large (infinite actually).\nThese observations lie right on the edge of the realistic range of the model for the data (maybe even over the edge).")

  }
  # if( type=="RQR.sim"){
  #   nsim <- 1000
  #   if( is.null( mc.cores))
  #     mc.cores <- getOption("mc.cores", 4)
  #   resids <- matrix( NA, nrow=object$n, ncol=object$S)
  #   RQR.fun <- function(ii){
  #     if( !quiet)
  #       setTxtProgressBar(pb, ii)
  #     X1 <- kronecker( matrix( 1, ncol=1, nrow=nsim), fm$titbits$X[ii,,drop=FALSE])
  #     W1 <- kronecker( matrix( 1, ncol=1, nrow=nsim), fm$titbits$W[ii,,drop=FALSE])
  #     sims <- regional_mix.simulate( nRCP=object$nRCP, S=object$S, n=nsim, p.x=object$p.x, p.w=object$p.w, alpha=object$coef$alpha, tau=object$coef$tau, beta=object$coef$beta, gamma=object$coef$gamma, logDisps=object$coef$disp, powers=object$titbits$power, X=X1, W=W1, offset=object$titbits$offset,family=object$family)
  #     sims <- sims[,1:object$S]
  #     yi <- object$titbits$Y[ii,,drop=FALSE]
  #     many_yi <- matrix( rep( yi, each=nsim), ncol=object$S)
  #     F_i <- colMeans( sims <= many_yi)
  #     F_i_minus <- colMeans( sims < many_yi)
  #     r_i <- runif( object$S, min=F_i_minus, max=F_i)
  #     return( qnorm( r_i))
  #   }
  #   if( !quiet)
  #     pb <- txtProgressBar(min = 1, max = object$n, style = 3, char = "><(('> ")
  #   if( Sys.info()['sysname'] == "Windows" | mc.cores==1)
  #     resids <- lapply( 1:object$n, RQR.fun)
  #   else
  #     resids <- parallel::mclapply( 1:object$n, RQR.fun, mc.cores=mc.cores)
  #   if( !quiet)
  #     message("")
  #   resids <- matrix( unlist( resids), nrow=object$n, ncol=object$S, byrow=TRUE)
  #   if( !quiet & sum( resids==Inf | resids==-Inf)>0)
  #     message( "Some residuals, well",sum( resids==Inf | resids==-Inf), "to be precise, are very large (infinite actually).\nThese observations lie right on the edge of the Monte Carlo approximation to the distribution function.\nThis may be remedied by getting a better approximation (increasing nsim).")
  # }
  return( resids)
}



#' @rdname regional_mix.simulate
#' @name regional_mix.simulate
#' @title Simulate a regional_mix dataset for modelling.
#' @description Simulates a data set from a mixture-of-experts model for RCP (for region of common profile) types.
#' @param nRCP Integer giving the number of RCPs
#' @param S Integer giving the number of species
#' @param n Integer giving the number of observations (sites)
#' @param p.x Integer giving the number of covariates (including the intercept) for the model for the latent RCP types
#' @param p.w Integer giving the number of covariates (excluding the intercept) for the model for the species data
#' @param alpha Numeric vector of length S. Specifies the mean prevalence for each species, on the logit scale
#' @param tau Numeric matrix of dimension c(nRCP-1,S). Specifies each species difference from the mean to each RCPs mean for the first nRCP-1 RCPs. The last RCP means are calculated using the sum-to-zero constraints
#' @param beta Numeric matrix of dimension c(nRCP-1,p.x). Specifies the RCP's dependence on the covariates (in X)
#' @param gamma Numeric matrix of dimension c(n,p.w). Specifies the species' dependence on the covariates (in W)
#' @param logDisps Logartihm of the (over-)dispersion parameters for each species for negative binomial, Tweedie and Normal models
#' @param powers Power parameters for each species for Tweedie model
#' @param X Numeric matrix of dimension c(n,p.x). Specifies the covariates for the RCP model. Must include the intercept, if one is wanted. Default is random numbers in a matrix of the right size.
#' @param W Numeric matrix of dimension c(n,p.w). Specifies the covariates for the species model. Must not include the intercept. Unless you want it included twice. Default is to give random levels of a two-level factor.
#' @param offset Numeric vector of size n. Specifies any offset to be included into the species level model.
#' @param family Text string. Specifies the family of the species data. Current options are "bernoulli" (default), "poisson", "negative.binomial", "tweedie" and "gaussian.
#' @export
#' @examples
#' \dontrun{
#' #generates synthetic data
#'set.seed( 151)
#'n <- 100
#'S <- 10
#'nRCP <- 3
#'my.dist <- "negative.binomial"
#'X <- as.data.frame( cbind( x1=runif( n, min=-10, max=10),
#'                           x2=runif( n, min=-10, max=10)))
#'Offy <- log( runif( n, min=30, max=60))
#'pols <- list()
#'pols[[1]] <- poly( X$x1, degree=3)
# Scale covariates so that regional_mix can get decent starting values
#'pols[[2]] <- poly( X$x2, degree=3)
#'X <- as.matrix( cbind( 1, X, pols[[1]], pols[[2]]))
#'colnames( X) <- c("const", 'x1', 'x2', paste( "x1",1:3,sep='.'),
#' paste( "x2",1:3,sep='.'))
#'p.x <- ncol( X[,-(2:3)])
#'p.w <- 3
#'W <- matrix(sample( c(0,1), size=(n*p.w), replace=TRUE), nrow=n, ncol=p.w)
#'colnames( W) <- paste( "w",1:3,sep=".")
#'alpha <- rnorm( S)
#'tau.var <- 0.5
#'b <- sqrt( tau.var/2)
#a double exponential for RCP effects
#'tau <- matrix( rexp( n=(nRCP-1)*S,
#' rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b), nrow=nRCP-1, ncol=S)
#'beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2, -76.0,
#' -14.2, -28.3, -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
#'gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
#'logDisp <- log( rexp( S, 1))
#'set.seed(121)
#'simDat <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n,
#' alpha=alpha, tau=tau, beta=beta, gamma=gamma, X=X[,-(2:3)], W=W,
#' family=my.dist, logDisp=logDisp, offset=Offy)
#'
#' }
"regional_mix.simulate" <- function(nRCP=3, S=20, n=200, p.x=3, p.w=0,
                                    alpha=NULL, tau=NULL, beta=NULL, gamma=NULL,
                                    logDisps=NULL, powers=NULL, X=NULL, W=NULL,
                                    offset=NULL, family="bernoulli"){
  if (is.null(alpha) | length(alpha) != S) {
    message("Random alpha from normal (-1,0.5) distribution")
    alpha <- rnorm(S,-1,0.5)
  }
  if (is.null(tau) | length(tau) != (nRCP - 1) * S) {
    message("Random tau from standard normal")
    tau <- rnorm( (nRCP-1)*S)
  }
  tau <- matrix(as.numeric(tau), nrow = nRCP - 1)
  if (is.null(beta) | length(beta) != (nRCP - 1) * p.x) {
    message("Random values for beta")
    beta <- rnorm( p.x*(nRCP-1))#as.numeric(c(0, 0, 0.4, 0, -0.2, 1))
  }
  beta <- matrix(as.numeric(beta), nrow = nRCP - 1)
  if( ( is.null(gamma) | length( gamma) != S * p.w)){
    if( p.w != 0){
      message("Random values for gamma")
      gamma <- rnorm( p.w*S)
      gamma <- matrix( as.numeric( gamma), nrow=S, ncol=p.w)
    }
    else
      gamma <- NULL
  }
  else
    gamma <- matrix( as.numeric( gamma), nrow=S)
  if( family == "negative.binomial" & (is.null( logDisps) | length( logDisps) != S)){
    message( "Random values for overdispersions")
    logDisps <- log( 1 + rgamma( n=S, shape=1, scale=0.75))
  }
  if( family=="gaussian" & (is.null( logDisps) | length( logDisps) != S)){
    message( "Random values for species' variance parameters")
    logDisps <- log( 1 + rgamma( n+S, shape=1, scale=0.75))
  }
  sppNames <- paste("spp", 1:S, sep = "")
  if (is.null(X)) {
    message("creating a RCP-level design matrix with random numbers")
    X <- cbind(1, matrix(runif(n * (p.x - 1), min = -10, max = 10), nrow = n))
    if( p.x > 1)
      colnames(X) <- c("intercept", paste("x", 1:(p.x - 1), sep = ""))
    else
      colnames(X) <- "intercept"
  }
  if( p.w>0)
    if( is.null( W)){
      message("Creating a species-level design matrix with random factor levels")
      W <- matrix(sample( c(0,1), size=(n*p.w), replace=TRUE), nrow=n, ncol=p.w)
      colnames(W) <- c(paste("w", 1:p.w, sep = ""))
    }
  if( is.null( offset))
    offset <- rep( 0, n)
  if( !family%in%c("bernoulli","poisson","negative.binomial","tweedie","gaussian")){
    message( "family not found, please choose from c('bernoulli','poisson','negative.binomial','tweedie','gaussian')")
    return( NA)
  }

  etaPi <- X %*% t(beta)
  pis <- t(apply(etaPi, 1, additive_logistic))
  habis <- apply(pis, 1, function(x) sample(1:nRCP, 1, FALSE, x))

  tau <- rbind(tau, -colSums(tau))
  etaMu <- tau + rep(alpha, each = nRCP)
  etaMu1 <- array( rep( offset, each=nRCP*S), dim=c(nRCP,S,n))
  if( p.w > 0){
    etaMu2 <- W %*% t( gamma)
    for( hh in 1:nRCP)
      etaMu1[hh,,] <- etaMu1[hh,,] + t( etaMu2)
  }
  for( hh in 1:nRCP)
    etaMu1[hh,,] <- etaMu1[hh,,] + rep( etaMu[hh,], times=n)
  etaMu <- etaMu1

  if( family=="bernoulli")
    mu <- inv.logit(etaMu)
  if( family %in% c("poisson","negative.binomial"))#,"tweedie"))
    mu <- exp( etaMu)
  if( family == "gaussian")
    mu <- etaMu

  fitted <- matrix( NA, nrow=n, ncol=S)
  for( ii in 1:n)
    fitted[ii,] <- mu[habis[ii], ,ii]

  if( family=="bernoulli")
    outcomes <- matrix(rbinom(n * S, 1, as.numeric( fitted)), nrow = n, ncol = S)
  if( family=="poisson")
    outcomes <- matrix(rpois(n * S, lambda=as.numeric( fitted)), nrow = n, ncol = S)
  if( family=="negative.binomial")
    outcomes <- matrix(rnbinom(n * S, mu=as.numeric( fitted), size=1/rep(exp( logDisps), each=n)), nrow = n, ncol = S)
  if( family=="gaussian")
    outcomes <- matrix( rnorm( n=n*S, mean=as.numeric( fitted), sd=rep( exp( logDisps), each=n)), nrow=n, ncol=S)

  colnames(outcomes) <- paste("spp", 1:S, sep = "")
  if( !all( offset==0))
    res <- as.data.frame(cbind(outcomes, X, W, offset))
  else
    res <- as.data.frame(cbind(outcomes, X, W))
  attr(res, "RCPs") <- habis
  attr(res, "pis") <- pis
  attr(res, "alpha") <- alpha
  attr(res, "tau") <- tau[-nRCP, ]
  attr(res, "beta") <- beta
  attr(res, "gamma") <- gamma
  attr(res, "logDisps") <- logDisps
  attr(res, "mu") <- mu
  return(res)
}

#'@title What is the average species profile per RCP?
#'@rdname regional_mix.species_profile
#'@name regional_mix.species_profile
#'@param object A RCP model
#'@param object2 A RCP model bootstrap object
#'@param CI The confidence intervals to report the range of values form bootstrap
#'@param \\dots Ignored for now.
#'@export
#'@description Extracts the average species' profile for each RCP.

"regional_mix.species_profile" <- function(object, object2=NULL, CI=c(0.025,0.975), type="response", ...){

  if(!type%in%c("response","link"))stop("prediction type not avaliable")
  if(is.null(object2)){
    if(check_if_sampling(object)) method <- "single_no_sp_results"
    else method <- "single_results"
  } else {
    method <- "bootstrap_results"
  }
  partial_mus <- switch(method,
                        single_no_sp_results = partial_mus_no_species_form(object,type),
                        single_results = partial_mus_with_species_form(object,type),
                        bootstrap_results = partial_mus_from_boostrap(object, object2, CI = CI,type))
  class(partial_mus) <- "regional_mix_membership"
  return(partial_mus)
}


#'@rdname stability.regional_mix
#'@name stability.regional_mix
#'@title Diagnostic checks to see if RCP groups are stable
#'@description For increasing size of hold-out samples, cooks distance and predictive log-likelihood are calculated and optionally plotted.
#'@param model a regional_mix model, as obtained by the function \code{regional_mix}. This is the model whose stability is assessed. Model must contain titbits (see ?regional_mix and particular attention to the argument titbits=TRUE)
#'@param oosSizeRange the size of the (successive) hold-out samples. If NULL (default), then a sequence of 10 sizes, from 1 to 0.2*model$n is used. The more numbers in this range, the slower the function will run.
#'@param times the number of hold-out samples to use. If times=model$n and oosSize is 1, then the sample contains each and every site. Otherwise, it is a sample of size times from the possible combinations of possible hold-out sets.
#'@param mc.cores the number of cores to farm the jobs out to.
#'@param quiet should the progress bar be displayed (bar for each oosSizeRange)
#'@param doPlot should the plots be produced? Default is that they should be.
#'@details The plots produced are: 1) leave-some-out Cook's distance (see \code{\link{cooks.distance.regional_mix}}) against holdout sample size; and 2) the predictive log-likelihood for times sites, against the holdout sample size.
#' In both plots, the values from the original model have been added to the plot.
#'@return Produces a regional_mix_stab object. This is a list with the oosSizeRnage, disty (the mean Cook's Distance for each subset size), nRCP, n, predlogls (log-likelihood of out-of-sample sites), logl.sites (the in-sample log-likelihood for full data set).
#'@export
#'@examples
#'\dontrun{
#'#not run as R CMD check complains about the time taken.
#'#This code will take a little while to run (about 3.5minutes on my computer)
#'  stability.regional_mix( fm, oosSizeRange=seq( from=1,to=fm$n,length=5),
#'                     times=fm$n, mc.cores=2, doPlot=FALSE);
#'}
#'
"stability.regional_mix" <- function( model, oosSizeRange=NULL, times=model$n, mc.cores=1, quiet=FALSE, doPlot=TRUE){
  if( is.null( oosSizeRange))
    oosSizeRange <- round( seq( from=1, to=model$n%/%5, length=10))
  if( any( oosSizeRange < 1))
    stop( "Silly number of RCPs. Specified range is: ", oosSizeRange, " and they should all be >= 1")
  disty <- matrix( NA, nrow=length( oosSizeRange), ncol=model$nRCP)
  predlogls <- array( NA, dim=c(length( oosSizeRange), model$n, times)) #matrix( NA, nrow=length( oosSizeRange), ncol=times)
  for( ii in oosSizeRange){
    tmp <- cooks.distance( model, oosSize=ii, times=times, mc.cores=mc.cores, quiet=quiet)
    disty[oosSizeRange==ii,] <- colMeans( abs( tmp$cooksD))
    predlogls[oosSizeRange==ii,,] <- tmp$predLogL
    #predlogls[oosSizeRange==ii,] <- colMeans( tmp$predLogL, na.rm=TRUE)
  }
  ret <- list( oosSizeRange=oosSizeRange, disty=disty, nRCP=model$nRCP,n=model$n, predlogls=predlogls, logl.sites=model$logl.sites)
  class( ret) <- "regional_mix_stab"

  if( doPlot)
    plot( ret)

  invisible( ret)
}

#' @rdname summary.regional_mix
#' @name summary.regional_mix
#' @title A summary from a regional_mix object.
#' @description A summary from a regional_mix object.
#' @param object A regional_mix model object
#' @param \\dots ignored.
#' @details A table is printed that contains the coefficient values, their standard errors, and their z-statistic. The second and thrid columns may be unreliable for some parameters.
#' @export

"summary.regional_mix" <-
  function (object, ...)
  {
    if (is.null(object$vcov)) {
      object$vcov <- matrix(NA, nrow = length(unlist(object$coef)),
                            ncol = length(unlist(object$coef)))
      stop("No variance matrix has been supplied")

    }
    message("Standard errors for alpha, tau and (probably) gamma parameters may be (are likely to be) misleading")
    res <- cbind(unlist(object$coefs), sqrt(diag(object$vcov)))
    res <- cbind(res, res[, 1]/res[, 2])
    res <- cbind(res, 2 * (1 - pnorm(abs(res[, 3]))))
    colnames(res) <- c("Estimate", "SE", "z-score", "p")
    return(res)
  }

#'@rdname vcov.regional_mix
#'@name vcov.regional_mix
#'@aliases vcov.regional_mix
#'@title Variance matrix for a regional_mix object.
#'@description Calculates variance-covariance matrix from a regional_mix object
#'@param object an object obtained from fitting a RCP (for region of common profile) mixture model. Such as that generated from a call to regional_mix(qv).
#'@param \\dots Other calls to the vcov function.
#'@param object2 an object of class \code{regional_mix} containing bootstrap samples of the parameter estimates (see regional_mix_boot(qv)). If NULL (default) the bootstrapping is performed from within the vcov function. If not null, then the vcov estimate is obtained from these bootstrap samples.
#'@param method the method to calculate the variance-covariance matrix. Options are:'FiniteDifference' (default), \code{BayesBoot}, \code{SimpleBoot}, and \code{EmpiricalInfo}. The two bootstrap methods (\code{BayesBoot} and \code{SimpleBoot}, see regional_mix_boot(qv)) should be more general and may possibly be more robust. The \code{EmpiricalInfo} method implements an empirical estimate of the Fisher information matrix, I can not recommend it however. It seems to behave poorly, even in well behaved simulations. It is computationally thrifty though.
#'@param nboot the number of bootstrap samples to take for the bootstrap estimation. Argument is ignored if !method \%in\% c(\code{FiniteDifference},'EmpiricalInfo').
#'@param mc.cores the number of cores to distrbute the calculations on. Default is 4. Set to 1 if the computer is running Windows (as it cannot handle forking -- see mclapply(qv)). Ignored if method=='EmpiricalInfo'.
#'@param D.accuracy The number of finite difference points used in the numerical approximation to the observed information matrix. Ignored if method != \code{FiniteDifference}. Options are 2 (default) and 4.
#'@details If method is \code{FiniteDifference}, then the estimates variance matrix is based on a finite difference approximation to the observed information matrix.
#'If method is either "BayesBoot" or "SimpleBoot", then the estimated variance matrix is calculated from bootstrap samples of the parameter estimates. See Foster et al (in prep) for details of how the bootstrapping is actually done, and regional_mix_boot(qv) for its implementation.
#'@return A square matrix of size equal to the number of parameters. It contains the variance matrix of the parameter estimates.
#'@references Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245
#' Foster, S.D., Hill, N.A. and Lyons, M., 2017. Ecological grouping of survey sites when sampling artefacts are present. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), pp.1031-1047.
#'@export
"vcov.regional_mix" <- function (object, ..., object2=NULL,
                                 method = "FiniteDifference",
                                 nboot = 1000, mc.cores=1, D.accuracy=2){

  if( method %in% c("simple","Richardson"))
    method <- "FiniteDifference"
  if (!method %in% c("FiniteDifference", "BayesBoot", "SimpleBoot", "EmpiricalInfo")) {
    error("Unknown method to calculate variance matrix, viable options are: 'FiniteDifference' (numerical), 'BayesBoot' (bootstrap), 'SimpleBoot' (bootstrap)', and 'EmpiricalInfo'.")
    return(NULL)
  }
  if( Sys.info()['sysname'] == "Windows")
    mc.cores <- 1
  X <- object$titbits$X
  p.x <- ncol( X)
  if( class( object$titbits$species_formula)=="formula"){
    form.W <- object$titbits$species_formula
    W <- object$titbits$W
    p.w <- ncol( W)
  }
  else{
    form.W <- NULL
    W <- -999999
    p.w <- 0
  }
  offy <- object$titbits$offset
  wts <- object$titbits$wts
  Y <- object$titbits$Y
  disty <- object$titbits$disty
  power <- object$titbits$power
  S <- object$S
  nRCP <- object$nRCP
  p.x <- object$p.x
  p.w <- object$p.w
  n <- object$n
  disty <- object$titbits$disty
  control <- object$titbits$control
  pis <- as.numeric( matrix( -999999, nrow = n, ncol = nRCP))
  mus <- as.numeric( array( -999999, dim=c( n, S, nRCP)))
  logCondDens <- as.numeric( matrix( -999999, nrow = n, ncol = nRCP))
  logls <- as.numeric(rep(-999999, n))
  alpha.score <- as.numeric(rep(-999999, S))
  tau.score <- as.numeric(matrix(-999999, nrow = nRCP - 1, ncol = S))
  beta.score <- as.numeric(matrix(-999999, nrow = nRCP - 1, ncol = p.x))
  if( p.w > 0)
    gamma.score <- as.numeric( matrix( -999999, nrow = S, ncol = p.w))
  else
    gamma.score <- -999999
  if( !is.null( object$coef$disp))
    disp.score <- as.numeric( rep( -999999, S))
  else
    disp.score <- -999999
  conv <- FALSE

  if (method %in% c("FiniteDifference")) {
    my.fun <- function(x) {
      start <- 0
      alpha <- x[start + 1:S]
      start <- start + S
      tau <- x[start + 1:((nRCP - 1) * S)]
      start <- start + (nRCP-1)*S
      beta <- x[start + 1:((nRCP - 1) * p.x)]
      start <- start + (nRCP-1)*p.x
      if( p.w > 0){
        gamma <- x[start + 1:(S*p.w)]
        start <- start + S*p.w
      }
      else
        gamma <- -999999
      if( any( !is.null( object$coef$disp)))
        disp <- x[start + 1:S]
      else
        disp <- -999999
      scoreContri <- -999999
      tmp <- .Call("RCP_C", as.numeric(Y), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
                   as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
                   alpha, tau, beta, gamma, disp, power,
                   as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
                   alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                   pis, mus, logCondDens, logls,
                   as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
                   as.integer( FALSE), as.integer( FALSE), as.integer( TRUE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")

      tmp1 <- c(alpha.score, tau.score, beta.score)
      if( p.w > 0)#class( object$titbits$species_formula) == "formula")
        tmp1 <- c( tmp1, gamma.score)
      if( !is.null( object$coef$disp))
        tmp1 <- c( tmp1, disp.score)
      return(tmp1)
    }
    hess <- nd2(x0=unlist( object$coefs), f=my.fun, mc.cores=mc.cores, D.accur=D.accuracy)#numDeriv::jacobian(my.fun, unlist(object$coefs), method = method)
    vcov.mat <- try( -solve(hess))
    if( inherits( vcov.mat, 'try-error')){
      attr(vcov.mat, "hess") <- hess
      warning( "Hessian appears to be singular and its inverse (the vcov matrix) cannot be calculated\nThe Hessian is returned as an attribute of the result (for diagnostics).\nMy deepest sympathies.  You could try changing the specification of the model, increasing the penalties, or getting more data.")
    }
    else
      vcov.mat <- ( vcov.mat + t(vcov.mat)) / 2 #to ensure symmetry
  }
  if( method %in% c( "BayesBoot","SimpleBoot")){
    object$titbits$control$optimise <- TRUE #just in case it was turned off (see regional_mix.multfit)
    if( is.null( object2))
      coefMat <- regional_mix_boot( object, nboot=nboot, type=method, mc.cores=mc.cores, quiet=TRUE)#, orderSamps=FALSE)
    else
      coefMat <- object2
    vcov.mat <- cov( coefMat)
  }
  if( method=="EmpiricalInfo"){
    message( "Information approximated by empirical methods.  I have not been able to get this to work, even for simulated data.  I hope that you are feeling brave!")
    alpha <- object$coef$alpha
    tau <- object$coef$tau
    beta <- object$coef$beta
    if( p.w > 0)
      gamma <- object$coef$gamma
    else
      gamma <- -999999
    if( any( !is.null( object$coef$disp)))
      disp <- object$coef$disp
    else
      disp <- -999999
    scoreContri <- as.numeric( matrix( NA, nrow=n, ncol=length( unlist( object$coef))))
    tmp <- .Call("RCP_C", as.numeric(Y), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
                 as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
                 alpha, tau, beta, gamma, disp, power,
                 as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
                 alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                 pis, mus, logCondDens, logls,
                 as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
                 as.integer( FALSE), as.integer( FALSE), as.integer( TRUE), as.integer( TRUE), as.integer( TRUE), PACKAGE = "ecomix")
    scoreContri <- matrix( scoreContri, nrow=n)
    summy <- matrix( 0, ncol=ncol( scoreContri), nrow=ncol( scoreContri))
    for( ii in 1:n){
      summy <- summy + scoreContri[ii,] %o% scoreContri[ii,]
    }
    tmp <- colSums( scoreContri)
    tmp <- tmp %o% tmp / n
    emp.info <- summy - tmp
    #    diag( emp.info) <- diag( emp.info) + 0.00001 #makes it invertable but not realistic.
    vcov.mat <- try( solve( emp.info))
    if( inherits( vcov.mat, 'try-error')){
      attr(vcov.mat, "hess") <- emp.info
      warning( "Empirical information matrix (average of the cross-products of the scores for each observation) appears to be singular and its inverse (the vcov matrix) cannot be calculated\nThe empirical inverse is returned as an attribute of the result (for diagnostics).\nMy deepest sympathies.  You could try changing the specification of the model, increasing the penalties, or getting more data. Note that you have chosen to use method=\"EmpricalInfo\", which is likely to cause heartache (albeit computationally thrifty heartache) -- try other methods (and probably do that first).")
    }
    else
      vcov.mat <- ( vcov.mat + t(vcov.mat)) / 2 #to ensure symmetry
  }

  return(vcov.mat)
}

