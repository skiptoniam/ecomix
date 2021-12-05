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
  if( object$p.w>0){
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

# "effect_data.regional_mix" <- function(focal.predictors, mod, ngrid = 50, ...){
#
#   if (is.null(mod$titbits))
#     stop("Model doesn't contain all information required for effectsPlotData.
#          Please supply model with titbits (from titbits=TRUE in species_mix call)")
#
#   Mode <- function(x, na.rm = FALSE) {
#     if (na.rm) {
#       x = x[!is.na(x)]
#     }
#     ux <- unique(x)
#     return(ux[which.max(tabulate(match(x, ux)))])
#   }
#
#   ## set up the data objects
#   X <- mod$titbits$X
#   W <- mod$titbits$W
#
#   ## set up the variables in the formula
#   tt <- terms(mod$titbits$rcp_formula)
#   tt <- delete.response(tt)
#   vars <- all.vars(parse(text=tt))
#   if(ncol(W)>1){
#     tt2 <- terms(mod$titbits$species_formula)
#     tt2 <- delete.response(tt2)
#     vars <- c(vars,all.vars(parse(text=tt2)))
#   }
#   nvars = length(vars)
#
#   pred.data <- mod$titbits$data
#   pred.data <- pred.data[,vars]
#
#   ## check for factors
#   factors <- NULL
#   for(ii in 1:nvars){
#     factors[ii] <- is.factor(pred.data[,vars[ii]]) | is.character(pred.data[,vars[ii]])
#   }
#
#   ## check for focal.predictors in pred.data
#   focal.ids <- lapply(focal.predictors, grep, colnames(pred.data))
#
#   # lists for data structures
#   mfs <- list() #catch model.frames
#   f.focal <- list()
#   v.focal <- list()
#   n.focal <- list()
#   for(i in 1:length(focal.ids)){
#
#     f.focal[[i]] <- factors[focal.ids[[i]]]
#     v.focal[[i]] <- pred.data[, focal.ids[[i]],drop=FALSE]
#     n.focal[[i]] <- seq_len(nvars)[-unlist(focal.ids[[i]])]
#
#     xx <- list()
#     for(j in 1:length(focal.ids[[i]])){
#       if(f.focal[[i]][j]) {
#         xx[[j]] = levels(v.focal[[i]][j])
#         ngrid = length(xx)
#       } else {
#         mi = min(v.focal[[i]][j])
#         ma = max(v.focal[[i]][j])
#         xx[[j]] = seq(mi, ma, length.out = ngrid)
#       }
#     }
#     XDataNew = data.frame(xx, stringsAsFactors = TRUE)
#     colnames(XDataNew) = vars[focal.ids[[i]]]
#     for (k in seq_len(length(n.focal[[i]]))) {
#       non.focal = n.focal[[i]][k]
#       f.non.focal = factors[non.focal]
#       v.non.focal = pred.data[, vars[non.focal]]
#       if (f.non.focal) {
#         XDataNew[, vars[non.focal]] = Mode(v.non.focal)
#       }
#       if (!f.non.focal) {
#         v.non.focal = pred.data[, vars[non.focal]]
#         XDataNew[, vars[non.focal]] = mean(v.non.focal)
#       }
#     }
#     mfs[[i]] <- XDataNew[,vars]
#   }
#
#   names(mfs) <- focal.predictors
#   class(mfs) <- "regional_mix_effectPlotData"
#
#   return(mfs)
#
# }


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

#' @rdname logLik.regional_mix
#' @name logLik.regional_mix
#' @title Extract log-likelihood from a regional_mix model.
#' @param object A fitted regional_mix model
#' @param \\dots Ignored
#' @export

"logLik.regional_mix" <-function (object, ...){
  return(object$logl)
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

# #' @rdname plot.regional_mix.sampling
# #' @name plot.regional_mix.sampling
# #' @title Plot residuals from a regional_mix model.
# #' @param object a fitted regional_mix model you wish to plot
# #' @param object2 a bootstrap object from regional_mix.bootstrap
# #' @param type What type of residuals to plot? 'RQR'; random quantile residuals or 'deviance' residuals.
# #' @param nsim Number of simulations to run
# #' @param alpha.conf The bounds of the confidence intervals.
# #' @param quiet Run in quiet mode.
# #' @param species Which species to plot as residuals.
# #' @param fitted.scale What scale to plot the residuals on?
# #' @details The two types of residuals are inherently different. The "RQR" residuals produce a residual for each species at each site and the "deviance" residuals produce a site residual (no species level residual). The plots also differ, the "RQR" type generates a single normal QQ-plot for all species and all sites, and a residual versus fitted plot for all species and sites (Described in Foster et al, 2013). The "deviance" type generates a pair of Tukey mean-difference plots, similar in spirit to a QQ-plot. The first is for point-wise confidence intervals and the second is for approximate global intervals. See Foster et al (2013) for details.
# #' The family for the "RQR" residuals should be standard normal. For "deviance" residuals, the distribution is unknown and simulation is used to graphically assess how odd the observed residuals look compared to ones generated assuming the model is correct.
# #' @references Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.
# #' Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245
# #' Foster, S.D., Hill, N.A. and Lyons, M., 2017. Ecological grouping of survey sites when sampling artefacts are present. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), pp.1031-1047.
# #' @export

# "plot.regional_mix.sampling" <-function(x, ... , object2=NULL, sampling_levels,
#                                         CI=c(0.025, 0.975), col="black", lty=1){
#
#   require(lattice)
#   gammas<-grepl("gamma",dimnames(object2)[[2]])
#   temp_dat<-object2[,gammas]
#
#   temp<-data.frame(avs=as.numeric(unname(colMeans(temp_dat))),
#                    t(apply(temp_dat, 2, quantile, probs=CI)),
#                    #sampling_var=rep(sampling_names, each=length(length(best_mod$names$spp))),
#                    sampling_var=sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 3),
#                    #Species=factor(rep(best_mod$names$spp,length(sampling_names))))
#                    Species=factor(sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 1)))
#
#   names(temp)[2:3]<-c("lower", "upper")
#   temp$Species<-gsub("."," ", temp$Species, fixed=TRUE) #get rid of '.' in species names
#   temp$Species<-as.factor(temp$Species) #convert back to factor
#   temp$Species <- factor(temp$Species, levels=rev(levels(temp$Species)))
#
#   trellis.par.set(superpose.symbol=list(pch=16,col=col, cex=1.2),
#                   superpose.line=list(col="transparent"))
#   dotplot(Species ~ avs, groups=sampling_var, data=temp, cols=col, lty=lty, low=temp$lower, high=temp$upper, subscript=TRUE,
#           auto.key=list(space="top", columns=2, cex=1.4, text=legend_fact),
#           ylab=list(cex=1.4), xlab=list("Coefficient",cex=1.4),
#           scales = list(tck = c(1, 0), x=list(cex=1.2), y=list(cex=1.2)),
#           prepanel = function(x, y, ...) { list(xlim=range(temp$lower, temp$upper)) },
#           panel=panel.superpose,
#           panel.groups=function(x, y, subscripts, group.number, cols, low, high, ...)
#           {
#             if(group.number==1) jiggle <- 0.1 else jiggle <- -0.1
#             panel.abline(v=0, lty=2)
#             panel.abline(h=1:length(best_mod$names$spp), col.line="light grey", lty=1)
#             panel.dotplot(x, y+jiggle, group.number, ...)
#             panel.arrows(low[subscripts], y+jiggle, high[subscripts], y+jiggle, code=3, angle=90,
#                          length=0.05, col=cols[group.number], lty=lty[group.number])
#             #panel.segments(temp$lower, y+jiggle, + temp$upper, y+jiggle, lty = lty, col =col, lwd=2, cex=1.2)
#           })
#
# }



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

#'
#
# "plot.species_effects" <-function(x,#best_mod,              # output of regimix function for final model
#                                 boot_obj,              # output of regiboot function for final model
#                                 legend_fact,           # levels of categorical sampling variable to plot. Coefficients relative to first level of factor.
#                                 # So usually 2:n levels(sampling_variable)
#                                CI=c(0.025, 0.975),    # confidence interval to plot
#                             col="black",           # colour/s of dots and CI lines. Specified in same way as col is usually specified
#                             lty=1)                 # lty= line type of CI lines. Specified in same way as lty is usually specified
# {
#
#   require(lattice)
#   # gammas<- paste0("gamma", 1: (length(Species)*length(sampling_names)))
#   gammas<-grepl("gamma",dimnames(boot_obj)[[2]])
#   temp_dat<-boot_obj[,gammas]
#
#   temp<-data.frame(avs=as.numeric(unname(colMeans(temp_dat))),
#                    t(apply(temp_dat, 2, quantile, probs=CI)),
#                    #sampling_var=rep(sampling_names, each=length(length(best_mod$names$spp))),
#                    sampling_var=sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 3),
#                    #Species=factor(rep(best_mod$names$spp,length(sampling_names))))
#                    Species=factor(sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 1)))
#
#   names(temp)[2:3]<-c("lower", "upper")
#   temp$Species<-gsub("."," ", temp$Species, fixed=TRUE) #get rid of '.' in species names
#   temp$Species<-as.factor(temp$Species) #convert back to factor
#   temp$Species <- factor(temp$Species, levels=rev(levels(temp$Species)))
#
#   trellis.par.set(superpose.symbol=list(pch=16,col=col, cex=1.2),
#                   superpose.line=list(col="transparent"))
#   dotplot(Species ~ avs, groups=sampling_var, data=temp, cols=col, lty=lty, low=temp$lower, high=temp$upper, subscript=TRUE,
#           auto.key=list(space="top", columns=2, cex=1.4, text=legend_fact),
#           ylab=list(cex=1.4), xlab=list("Coefficient",cex=1.4),
#           scales = list(tck = c(1, 0), x=list(cex=1.2), y=list(cex=1.2)),
#           prepanel = function(x, y, ...) { list(xlim=range(temp$lower, temp$upper)) },
#           panel=panel.superpose,
#           panel.groups=function(x, y, subscripts, group.number, cols, low, high, ...)
#           {
#             if(group.number==1) jiggle <- 0.1 else jiggle <- -0.1
#             panel.abline(v=0, lty=2)
#             panel.abline(h=1:length(best_mod$names$spp), col.line="light grey", lty=1)
#             panel.dotplot(x, y+jiggle, group.number, ...)
#             panel.arrows(low[subscripts], y+jiggle, high[subscripts], y+jiggle, code=3, angle=90,
#                          length=0.05, col=cols[group.number], lty=lty[group.number])
#             #panel.segments(temp$lower, y+jiggle, + temp$upper, y+jiggle, lty = lty, col =col, lwd=2, cex=1.2)
#           })
#
# }

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
#' set.seed( 151)
#' n <- 100
#' S <- 10
#' nRCP <- 3
#' my.dist <- "negative.binomial"
#' X <- as.data.frame( cbind( x1=runif( n, min=-10, max=10),
#'                           x2=runif( n, min=-10, max=10)))
#' Offy <- log( runif( n, min=30, max=60))
#' pols <- list()
#' pols[[1]] <- poly( X$x1, degree=3)
#' pols[[2]] <- poly( X$x2, degree=3)
#' X <- as.matrix( cbind( 1, X, pols[[1]], pols[[2]]))
#' colnames( X) <- c("const", 'x1', 'x2', paste( "x1",1:3,sep='.'),
#' paste( "x2",1:3,sep='.'))
#' p.x <- ncol( X[,-(2:3)])
#' p.w <- 3
#' W <- matrix(sample( c(0,1), size=(n*p.w), replace=TRUE), nrow=n, ncol=p.w)
#' colnames( W) <- paste( "w",1:3,sep=".")
#' alpha <- rnorm( S)
#' tau.var <- 0.5
#' b <- sqrt( tau.var/2)
#' tau <- matrix( rexp( n=(nRCP-1)*S,rate=1/b) - rexp( n=(nRCP-1)*S, rate=1/b),
#'  nrow=nRCP-1, ncol=S)
#' beta <- 0.2 * matrix( c(-1.2, -2.6, 0.2, -23.4, -16.7, -18.7, -59.2,
#'  -76.0,-14.2, -28.3, -36.8, -17.8, -92.9,-2.7), nrow=nRCP-1, ncol=p.x)
#' gamma <- matrix( rnorm( S*p.w), ncol=p.w, nrow=S)
#' logDisp <- log( rexp( S, 1))
#' set.seed(121)
#' simDat <- regional_mix.simulate( nRCP=nRCP, S=S, p.x=p.x, p.w=p.w, n=n,
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
#'@param type What type of prediction to perform? Default is response, alternative is 'link'.
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
  class(partial_mus) <- "regional_mix_profile"
  return(partial_mus)
}

"partial_mus_no_species_form" <- function(object, type, ...){

  ## what are the species taus?
  tau <- coef(object)$tau
  tau <- rbind(tau, -colSums( tau))

  ## what was the the model offset?
  offy <- object$titbits$offset

  ## what is the linear predictor (eta)
  eta <- sweep(tau, 2, object$coefs$alpha, "+") + mean(offy)

  ## what is the link function of appropriate family?
  if(object$family=="bernoulli")
    link.fun <- stats::make.link('logit')
  if(object$family%in%c("poisson","negative.binomial"))
    link.fun <- stats::make.link('log')
  if(object$family=='guassian')
    link.fun <- stats::make.link('identity')

  ## what are the partial mus: dim[nRCPs,nSpp]
  if(type%in%"response")partial_mus <- link.fun$linkinv(eta)
  else partial_mus <- eta
  dimnames(partial_mus)[[2]] <- object$names$spp
  dimnames(partial_mus)[[1]] <- object$names$RCPs

  ## return the partial mus if their is no sampling artifacts (species formula).
  return(partial_mus)
}


"partial_mus_with_species_form" <- function(object, type, ... ){

  ## what are the taus?
  tau <- coef(object)$tau
  tau <- rbind(tau, -colSums( tau))

  ## what was the the model offset?
  offy <- object$titbits$offset

  ## what is the linear predictor (eta)?
  eta <- sweep(tau, 2, coef(object)$alpha, "+") + mean(offy)

  ## what is the link function of appropriate family?
  if(object$family=="bernoulli")
    link.fun <- stats::make.link('logit')
  if(object$family%in%c("poisson","negative.binomial"))
    link.fun <- stats::make.link('log')
  if(object$family=='guassian')
    link.fun <- stats::make.link('identity')

  ## probabilities for all other levels of same sampling var
  res<- lapply(seq_along(object$names$Wvars),function(jj){
    new_eta <- sweep(eta, 2, coef(object)$gamma[,object$names$Wvars[jj]], "+");
    if(type%in%"response")part_mus <- link.fun$linkinv(eta)
    else part_mus <- eta
    return(part_mus)})

  names(res)<- object$names$Wvars
  return(res)
}

"partial_mus_from_boostrap"  <- function(object, object2, CI=c(0.025,0.975), type){

  #set up coefficient extraction
  taus<-grepl("tau",dimnames(object2)[[2]])
  alphas<-grepl("alpha",dimnames(object2)[[2]])

  ## what is the link function of appropriate family?
  if(object$family=="bernoulli") link.fun <- stats::make.link('logit')
  if(object$family%in%c("poisson","negative.binomial")) link.fun <- stats::make.link('log')
  if(object$family=='negative.binomial') link.fun <- stats::make.link('log')
  if(object$family=='guassian') link.fun <- stats::make.link('identity')

  if(check_if_sampling(object)){

    res_all <- list()
    for(i in seq_len(dim(object2)[1])){

      ## bootstrap alpha (intercept)
      tmp_alphas<-object2[i,alphas]

      # bootstrap tau
      tmp_tau <- object2[i,taus]
      tmp_tau <- matrix(tmp_tau, nrow=length(object$names$RCPs)-1)
      tmp_tau_all <- rbind(tmp_tau,-colSums(tmp_tau))
      colnames(tmp_tau_all) <- object$names$spp
      rownames(tmp_tau_all) <- object$names$RCPs

      ## offset from the model if used.
      offy <- object$titbits$offset

      ## what is the linear predictor (eta)
      tmp_eta <- sweep(tmp_tau_all, 2, tmp_alphas, "+") + mean(offy)

      #calculate values
      if(type%in%"response")part_mus <- link.fun$linkinv(tmp_eta)
      if(type%in%"link") part_mus <- tmp_eta
      res_all[[i]]<-as.matrix(part_mus)
    }

    overall_temp<-array(unlist(res_all), dim=c( length(object$names$RCPs),length(object$names$spp),nrow(object2)))
    overall_res<-list( mean=apply(overall_temp, c(1,2), mean),
                       sd= apply(overall_temp, c(1,2), sd),
                       lower= apply(overall_temp, c(1,2), function(x) quantile(x, probs=CI[1])),
                       upper= apply(overall_temp, c(1,2), function(x) quantile(x, probs=CI[2])))

    dimnames(overall_res[[1]])<-dimnames(overall_res[[2]])<-dimnames(overall_res[[3]])<-dimnames(overall_res[[4]])<-list(object$names$RCPs, object$names$spp)
    return (overall_res)
  }

  if (!check_if_sampling(object)){

    gammas<-grepl("gamma",dimnames(object2)[[2]])
    res_all <- list()
    # res <- rep( list(list()), length(object$names$Wvars))


    for(i in seq_len(dim(object2)[1])){
      ## bootstrap alpha (intercept)
      tmp_alphas<-object2[i,alphas]

      # bootstrap tau
      tmp_tau <- object2[i,taus]
      tmp_tau <- matrix(tmp_tau, nrow=length(object$names$RCPs)-1)
      tmp_tau_all <- rbind(tmp_tau,-colSums(tmp_tau))
      colnames(tmp_tau_all) <- object$names$spp
      rownames(tmp_tau_all) <- object$names$RCPs

      ## offset from the model if used.
      offy <- object$titbits$offset

      #gamma
      tmp_gamma<-object2[i, gammas]
      tmp_gamma<-matrix(tmp_gamma, nrow=length(object$names$spp))
      colnames(tmp_gamma)<-object$names$Wvars
      rownames(tmp_gamma)<-object$names$spp

      ## what is the linear predictor (eta)
      tmp_eta <- sweep(tmp_tau_all, 2, tmp_alphas, "+") + mean(offy)

      res<- lapply(seq_along(object$names$Wvars),function(jj){
        new_eta <- sweep(tmp_eta, 2, tmp_gamma[,jj], "+");
        if(type%in%"response")part_mus <- link.fun$linkinv(tmp_eta)
        if(type%in%"link") part_mus <- tmp_eta
        return(part_mus)})

      names(res)<-object$names$Wvars
      res_all[[i]] <- res
    }

    #Compile list of summaries at the sampling factor level
    samp_res <- rep(list(list()), length(object$names$Wvars))
    names(samp_res) <- object$names$Wvars

    for(k in seq_along(object$names$Wvars)){
      samp_res[[k]]<-list(mean=apply(simplify2array(res_all[[k]]), c(1,2), mean),
                          sd=apply(simplify2array(res_all[[k]]), c(1,2), sd),
                          lower=apply(simplify2array(res_all[[k]]), c(1,2), function(x) quantile(x, probs=CI[1])),
                          upper=apply(simplify2array(res_all[[k]]), c(1,2), function(x) quantile(x, probs=CI[2])))
    }

    overall_temp<-list()
    for(i in seq_len(dim(object2)[1])){
      get_vals <- res_all[[i]]
      overall_temp[[i]]<-apply(simplify2array(get_vals), c(1,2), mean)
    }

    overall_samp <-list(mean=apply(simplify2array(overall_temp), c(1,2), mean),
                        sd= apply(simplify2array(overall_temp), c(1,2), sd),
                        lower= apply(simplify2array(overall_temp), c(1,2), function(x) quantile(x, probs=CI[1])),
                        upper= apply(simplify2array(overall_temp), c(1,2), function(x) quantile(x, probs=CI[2])))
    samp_res$overall<-overall_samp
    return(samp_res)
  }
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
#'@param object2 an object of class \code{regional_mix} containing bootstrap samples of the parameter estimates (see regional_mix.bootstrap(qv)). If NULL (default) the bootstrapping is performed from within the vcov function. If not null, then the vcov estimate is obtained from these bootstrap samples.
#'@param method the method to calculate the variance-covariance matrix. Options are:'FiniteDifference' (default), \code{BayesBoot}, \code{SimpleBoot}, and \code{EmpiricalInfo}. The two bootstrap methods (\code{BayesBoot} and \code{SimpleBoot}, see regional_mix.bootstrap(qv)) should be more general and may possibly be more robust. The \code{EmpiricalInfo} method implements an empirical estimate of the Fisher information matrix, I can not recommend it however. It seems to behave poorly, even in well behaved simulations. It is computationally thrifty though.
#'@param nboot the number of bootstrap samples to take for the bootstrap estimation. Argument is ignored if !method \%in\% c(\code{FiniteDifference},'EmpiricalInfo').
#'@param mc.cores the number of cores to distrbute the calculations on. Default is 4. Set to 1 if the computer is running Windows (as it cannot handle forking -- see mclapply(qv)). Ignored if method=='EmpiricalInfo'.
#'@param D.accuracy The number of finite difference points used in the numerical approximation to the observed information matrix. Ignored if method != \code{FiniteDifference}. Options are 2 (default) and 4.
#'@details If method is \code{FiniteDifference}, then the estimates variance matrix is based on a finite difference approximation to the observed information matrix.
#'If method is either "BayesBoot" or "SimpleBoot", then the estimated variance matrix is calculated from bootstrap samples of the parameter estimates. See Foster et al (in prep) for details of how the bootstrapping is actually done, and regional_mix.bootstrap(qv) for its implementation.
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
      coefMat <- bootstrap( object, nboot=nboot, type=method, mc.cores=mc.cores, quiet=TRUE)#, orderSamps=FALSE)
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

