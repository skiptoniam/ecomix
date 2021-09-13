##### S3 class SAM functions #####
#' @rdname AIC.species_mix
#' @name AIC.species_mix
#' @title Return AIC from a species_mix model
#' @param object A species mix object
#' @param k AIC parameter
#' @param \\dots Ignored
#' @export

"AIC.species_mix" <- function (object, k=NULL, ...){
  effect.param <- length(object$beta) + object$G + object$S + object$npw*object$S +
    object$npu + ifelse(object$family %in% c("negative.binomial","tweedie","gaussian"),object$S,0)
  p <- effect.param
  if (is.null(k))
    k <- 2
  star.ic <- -2 * object$logl + k * p
  return(star.ic)
}

#' @rdname BIC.species_mix
#' @name BIC.species_mix
#' @title Return BIC from a species_mix model
#' @param object A species mix object
#' @param \\dots Ignored
#' @export

"BIC.species_mix" <-  function (object, ...){
  effect.param <- length(object$beta) + object$G + object$S + object$npw*object$S +
    object$npu + ifelse(object$family %in% c("negative.binomial","tweedie","gaussian"),object$S,0)
  p <- effect.param
  k <- log(object$n)
  star.ic <- -2 * object$logl + k * p
  return(star.ic)
}

#' @rdname coef.species_mix
#' @name coef.species_mix
#' @title print the coefficients from the species_mix model
#' @param object A species mix object
#' @param \\dots Ignored
#' @export
"coef.species_mix" <- function (object, ...){
  res <- list()
  res$alpha <- object$coefs$alpha
  names(res$alpha) <- object$names$spp
  if( !is.null( object$coef$beta)){
    res$beta <- matrix(object$coefs$beta, nrow = object$G, ncol = object$npx)
    rownames(res$beta) <- object$names$SAMs
    colnames(res$beta) <- object$names$Xvars
  }
  if((object$npw)>0){
    res$gamma <- matrix(object$coefs$gamma, nrow = object$S, ncol = object$npw)
    rownames(res$gamma) <- object$names$spp
    colnames(res$gamma) <- object$names$Wvars
  }
  if((object$npu)>0){
    res$delta <- object$coefs$delta
    names(res$delta) <- object$names$Uvars
  }
  if(object$family%in%c('negative.binomial','tweedie','gaussian')){
    res$theta <- object$coef$theta
    names(res$theta) <- object$names$spp
  }
  return(res)
}


#' @rdname logLik.species_mix
#' @name logLik.species_mix
#' @title Extract log-likelihood from a species_mix model.
#' @param object A fitted species_mix model
#' @param \\dots Ignored
#' @export

"logLik.species_mix" <-function (object, ...){
  return(object$logl)
}

#' @rdname plot.species_mix
#' @name plot.species_mix
#' @title plot.species_mix
#' @param x a fitted species_mix model.
#' @param species which species residuals to plot. Default is "AllSpecies".
#' @param type "response" or "link" for plotting of residuals to be on the linear predictor scale.
#' @param \\dots Extra plotting arguments.
#' @details Plot random quantile residuals (RQR). "RQR" produces residuals for each species.
#' @references  Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.
#' @export

"plot.species_mix" <- function (x,
                                species="AllSpecies",
                                type="response",
                                ...){
  # if( ! type %in% c("RQR"))
  # stop( "Unknown type of residuals. Options are 'RQR'.\n")
  if( ! all( species %in% c("AllSpecies",x$names$spp)))
    stop( "Unknown species.  Options are 'AllSpecies' or any one of the species names as supplied (and stored in x$names$spp)")

  # if( type=="RQR"){
  obs.resid <- residuals(x, type="RQR", quiet=TRUE)
  S <- x$S
  sppID <- rep( TRUE, S)
  if( species != "AllSpecies"){
    sppID <- x$names$spp %in% species
    obs.resid <- obs.resid[,sppID, drop=FALSE]
    S <- ncol( obs.resid)
  }
  if( sum( obs.resid==Inf | obs.resid==-Inf) > 0){
    message( "Infinite residuals removed from residual plots: ", sum( obs.resid==Inf | obs.resid==-Inf), " in total.")
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
  abline( 0,1,lwd=2)
  preds <- sam_internal_pred_species(x$coef$alpha, x$coef$beta, x$tau,
                                     x$coef$gamma, x$coef$delta, x$G, x$S, x$titbits$X,
                                     x$titbits$W,x$titbits$U, x$titbits$offset, x$family,
                                     type=type)

  preds <- preds[,sppID]

  # switch( fitted.scale,
  #         log = { loggy <- "x"},
  #         logit = { loggy <- ""; preds <- log( preds / (1-preds))},
  #         {loggy <- ""})
  plot( preds, obs.resid, xlab="Fitted", ylab="RQR", main="Residual versus Fitted", sub="Colours separate species", pch=20, col=rep( 1:S, each=x$n))
  abline( h=0)

  # }
}

#' @rdname plot.species_mix.multifit
#' @name plot.species_mix.multifit
#' @title plot.species_mix.multifit
#' @param x a fitted species_mix.multifit object.
#' @param type What type of response to plot options are "BIC", "AIC" or "logLik".
#' @param \\dots Extra plotting arguments.
#' @details Plot BIC or other likelihood based indices from multiple fit object.
#' @export
#' @importFrom stats BIC AIC

"plot.species_mix.multifit" <- function(x, type=c("BIC","AIC","logLik"), ...){

  if(x$groupselection){
  grps <- x$nArchetypes
  nSAMs_fm <- x$multiple_fits
  SAMsamp_ll <- sapply( nSAMs_fm, function(x) sapply( x, function(y) y$logl))
  SAMsamp_AICs <- sapply( nSAMs_fm, function(x) sapply( x, function(y) AIC(y)))
  SAMsamp_BICs <- sapply( nSAMs_fm, function(x) sapply( x, function(y) BIC(y)))
  SAMsamp_Gs <- sapply( nSAMs_fm, function(x) sapply( x, function(y) y$G))
  SAMsamp_minPosteriorSites <- sapply( nSAMs_fm, function(y) sapply( y, function(x) min( colSums( x$tau))))
  SAMsamp_ObviouslyBad <- which(SAMsamp_minPosteriorSites < 0.5)

  #ll
  SAMsamp_ll[SAMsamp_ObviouslyBad] <- NA
  SAMsamp_ll <- matrix(SAMsamp_ll, nrow = length(nSAMs_fm[[1]]))
  SAMsamp_minll <- apply( SAMsamp_ll, 2, min, na.rm=TRUE)

  # AIC
  SAMsamp_AICs[SAMsamp_ObviouslyBad] <- NA
  SAMsamp_AICs <- matrix(SAMsamp_AICs, nrow = length(nSAMs_fm[[1]]))
  SAMsamp_minAICs <- apply( SAMsamp_AICs, 2, min, na.rm=TRUE)

  #BIC
  SAMsamp_BICs[SAMsamp_ObviouslyBad] <- NA
  SAMsamp_BICs <- matrix(SAMsamp_BICs, nrow = length(nSAMs_fm[[1]]))
  SAMsamp_minBICs <- apply( SAMsamp_BICs, 2, min, na.rm=TRUE)

  type <- match.arg(type)

  if(type=="AIC"){
    df2a <- data.frame(grps=grps,indice=SAMsamp_minAICs)
    df2b <- data.frame(grps=rep( grps, each=nrow( SAMsamp_AICs)),indice=as.numeric(SAMsamp_AICs))
  }
  if(type=="BIC"){
  df2a <- data.frame(grps=grps,indice=SAMsamp_minBICs)
  df2b <- data.frame(grps=rep( grps, each=nrow( SAMsamp_BICs)),indice=as.numeric(SAMsamp_BICs))
  }
  if(type=="logLik"){
    df2a <- data.frame(grps=grps,indice=SAMsamp_minll)
    df2b <- data.frame(grps=rep( grps, each=nrow( SAMsamp_ll)),indice=as.numeric(SAMsamp_ll))
  }

  plot(indice~grps,
       data = df2a,
       pch=16,
       ylab=type,
       xlab="nArchetypes",
       ...)
  lines(indice~grps,
        data = df2a)
  points(indice~grps,data=df2b,pch=16)

  } else {
    message('Cannot plot because you have not done group selection.')
  }

}


#'@rdname predict.species_mix
#'@name predict.species_mix
#'@title Predict a species_mix model.
#'@param object is a matrix model returned from the species_mix model.
#'@param object2 is a species mix bootstrap object.
#'@param newdata a matrix of new observations for prediction.
#'@param offset an offset for prediction
#'@param nboot Number of bootstraps (or simulations if using IPPM) to run if no object2 is provided.
#'@param alpha confidence level. default is 0.95
#'@param mc.cores number of cores to use in prediction. default is 1.
#'@param type Do you want to predict the 'response' or the 'link'; ala glm style predictions.
#'@param prediction.type Do you want to produce 'archetype' or 'species' level predictions. default is 'archetype'.
#'@param na.action The type of action to apply to NA data. Default is "na.pass" see predict.lm for more details.
#'@param \\dots Ignored
#'@description Predict species archetypes from a species_mix model. You can also predict the conditional species predictions using "prediction.type='species'".
#'@export
#'@examples
#'\donttest{
#' library(ecomix)
#' set.seed(42)
#' sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),
#' collapse = ','),")~x1+x2"))
#' sp_form <- ~ 1
#' beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,1.9),3,2,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),
#' x2=stats::rnorm(100,0,2.5))
#' dat[,-1] <- scale(dat[,-1])
#' simulated_data <- species_mix.simulate(archetype_formula = sam_form,species_formula = sp_form,
#' data = dat,beta=beta,family="bernoulli")
#' fm1 <- species_mix(archetype_formula = sam_form,species_formula = sp_form,
#' data = simulated_data, family = 'bernoulli',  nArchetypes=3)
#' preds_fm1 <- predict(fm1)
#'}

"predict.species_mix" <- function(object, object2 = NULL, newdata = NULL,
                                  offset = NULL, nboot = 0, alpha = 0.95,
                                  mc.cores = 1, type = 'response',
                                  prediction.type='archetype',
                                  na.action = "na.pass", ...){
  if (is.null(newdata)) {
    X <- object$titbits$X
    W <- object$titbits$W
    U <- object$titbits$U
    offset <- object$titbits$offset
  } else {

    ## terms
    tt <- terms(object)

    ## Set up X based on arch.terms
    arch.tm <- tt[[1]]
    dat.levels <- lapply(newdata,levels)
    mfx <- model.frame(arch.tm, newdata, xlev = dat.levels)
    if (!is.null(cl <- attr(arch.tm, "dataClasses")))
      .checkMFClasses(cl, mfx)
    dat.fac <- vapply(newdata, is.factor, logical(1L))
    contrasts.list <- lapply(dat.fac,function(x)ifelse(x==TRUE,"contr.treatment",NA))
    contrasts.list <- Filter(Negate(anyNA),contrasts.list)
    X <- model.matrix(arch.tm, mfx, contrasts.arg = contrasts.list)
    X <- delete.intercept(X)

    ## Setup species matrix based on spp.terms
    spp.tm <- tt[[2]]
    if(length(attr(spp.tm,"factors"))>0){
      dat.levels <- lapply(newdata,levels)
    } else {
      dat.levels <- NULL
    }
    mfw <- model.frame(spp.tm, newdata, xlev = dat.levels)
    if (!is.null(cl <- attr(spp.tm, "dataClasses")))
      .checkMFClasses(cl, mfw)
    dat.fac <- vapply(newdata, is.factor, logical(1L))
    contrasts.list <- lapply(dat.fac,function(x)ifelse(x==TRUE,"contr.treatment",NA))
    contrasts.list <- Filter(Negate(anyNA),contrasts.list)
    W <- model.matrix(spp.tm, mfw, contrasts.arg = contrasts.list)

    if(!is.null(object$titbits$U)){
      all.tm <- tt[[3]]
      if(length(attr(all.tm,"factors"))>0){
        dat.levels <- lapply(newdata,levels)
      } else {
        dat.levels <- NULL
      }
      mfu <- model.frame(all.tm, newdata, xlev = dat.levels)
      if (!is.null(cl <- attr(all.tm, "dataClasses")))
        .checkMFClasses(cl, mfu)
      dat.fac <- vapply(newdata, is.factor, logical(1L))
      contrasts.list <- lapply(dat.fac,function(x)ifelse(x==TRUE,"contr.treatment",NA))
      contrasts.list <- Filter(Negate(anyNA),contrasts.list)
      U <- model.matrix(all.tm, mfu, contrasts.arg = contrasts.list)
     } else {
      U <- NULL
    }

    offset <- model.frame(arch.tm, data = newdata)
    offset <- model.offset(offset)
  }

  if (is.null(offset))
    offset <- rep(0, nrow(X))

  S <- object$S
  G <- object$G
  n <- object$n
  npx <- object$npx
  npw <- object$npw
  npu <- object$npu

  spp_wts <- object$titbits$spp_weights
  site_spp_wts <- object$titbits$site_spp_weights

  disty_cases <- c("bernoulli","poisson","ippm","negative.binomial","tweedie","gaussian","binomial")
  disty <- get_family_sam(disty_cases, object$titbits$family)
  tau <- object$tau
  if (is.null(object2)) {
    if (nboot > 0) {
      if( !object$titbits$control$quiet)
        message("Using a parametric bootstrap based on the ML estimates and their vcov")
      my.nboot <- nboot
    }
    else
      my.nboot <- 0
    allCoBoot <- species_mix_boot_parametric(object = object, nboot = my.nboot)
  } else {
    if( !object$titbits$control$quiet)
      message("Using supplied species_mix.bootstrap object (non-parametric bootstrap)")
    allCoBoot <- as.matrix(object2)
    nboot <- nrow(object2)
  }
  if (is.null(allCoBoot))
    return(NULL)

  alphaBoot <- allCoBoot[, seq_len(S), drop=FALSE]
  betaBoot <- allCoBoot[, S + seq_len((G*npx)), drop=FALSE]
  gammaBoot <- allCoBoot[, S + (G-1) + (G*npx) + seq_len((S*npw)), drop=FALSE]
  deltaBoot <- allCoBoot[, S + (G-1) + (G*npx) + (S*npw) + seq_len(npu), drop=FALSE]


  alphaIn <- c(NA, as.numeric(object$coefs$alpha))
  alphaIn <- alphaIn[-1]
  betaIn <- c(NA, as.numeric(object$coef$beta))
  betaIn <- betaIn[-1]
  # etaIn <- c(NA, as.numeric(object$coef$eta))
  # etaIn <- etaIn[-1]
  if (npw>0) {
    gammaIn <- c(NA, as.numeric(object$coef$gamma))
    gammaIn <- gammaIn[-1]
    usegamma <- 1
  } else {
    gammaIn <- -999999
    usegamma <- 0
  }
  if (npw>0) {
    deltaIn <- c(NA, as.numeric(object$coef$gamma))
    deltaIn <- deltaIn[-1]
    usedelta <- 1
  } else {
    deltaIn <- -999999
    usedelta <- 0
  }
  if (disty%in%c(4,6)) {
    thetaIn <- c(NA, as.numeric(object$coef$theta))
    thetaIn <- thetaIn[-1]
    usetheta <- 1
  } else {
    thetaIn <- -999999
    usetheta <- 0
  }

  outcomes <- matrix(NA, nrow = nrow(X), ncol = S)
  myContr <- object$titbits$control
  nam <- paste("G", 1:G, sep = "_")

  boot.funny.sam <- function(seg) {
    if (any(segments <= 1)) {
      nboot <- 0
      bootSampsToUse <- 1
      tmp <- switch (prediction.type,
                     archetype = sam_internal_pred_groups(alpha = object$coefs$alpha,
                                                          beta = object$coefs$beta,
                                                          gamma = object$coefs$gamma,
                                                          delta = object$coefs$delta,
                                                          tau = tau, G = G, S = S, X = X, W = W, U = U,
                                                          offset = offset, family = object$family, type = type),
                     species = sam_internal_pred_species(alpha = object$coefs$alpha,
                                                         beta = object$coefs$beta,
                                                         gamma = object$coefs$gamma,
                                                         delta = object$coefs$delta,
                                                         tau = tau, G = G, S = S, X = X, W = W,  U = U,
                                                         offset = offset, family = object$family, type = type))
    } else {
      nboot <- segments[seg]
      bootSampsToUse <- (sum( segments[1:seg])-segments[seg]+1):sum(segments[1:seg])

      # add in species level preds.
      tmp <- lapply(bootSampsToUse,function(ii)switch (prediction.type,
                                                       archetype = sam_internal_pred_groups(alpha = alphaBoot[ii,],
                                                                                            beta = matrix(betaBoot[ii,],G,npx),
                                                                                            gamma = matrix(gammaBoot[ii,],S,npw),
                                                                                            delta = deltaBoot[ii,],
                                                                                            tau = tau, G = G, S = S, X = X, W = W, U=U,
                                                                                            offset = offset, family = object$family,
                                                                                            type = type),
                                                       species = sam_internal_pred_species(alpha = alphaBoot[ii,],
                                                                                           beta = matrix(betaBoot[ii,],G,npx),
                                                                                           gamma = matrix(gammaBoot[ii,],S,npw),
                                                                                           delta = deltaBoot[ii,],
                                                                                           tau = tau, G = G, S = S, X = X, W = W, U=U,
                                                                                           offset = offset, family = object$family,
                                                                                           type = type)))

    }

    if (nboot == 0) {
      ret_grp <- tmp
      if(prediction.type%in%"archetype")colnames(ret_grp) <- object$names$SAMs
      if(prediction.type%in%"species")colnames(ret_grp) <- object$names$spp
      return(ret_grp)
    }

    if(prediction.type%in%"archetype"){
      bootPreds <- matrix(do.call("cbind",lapply(tmp,c)), nrow = nrow(X) * G,  ncol = nboot)
    }
    if(prediction.type%in%"species"){
      bootPreds <- matrix(do.call("cbind",lapply(tmp,c)), nrow = nrow(X) * S,  ncol = nboot)
    }
    return(bootPreds)
  }

  segments <- -999999
  ret <- list()
  ptPreds <- boot.funny.sam(1)
  if (nboot > 0) {
    if (Sys.info()["sysname"] == "Windows") {
      if( !object$titbits$control$quiet)
        message("Parallelised version of function not available for Windows machines. Reverting to single processor.")
      mc.cores <- 1
    }
    segments <- rep(nboot%/%mc.cores, mc.cores)
    if( nboot %% mc.cores > 0)
      segments[1:(nboot%%mc.cores)] <- segments[1:(nboot%%mc.cores)] + 1

    tmp <- parallel::mclapply(1:mc.cores, boot.funny.sam, mc.cores = mc.cores)
    bootPreds <- do.call("cbind", tmp)
    bPreds <- list()
    row.exp <- rowMeans(bootPreds)
    if(prediction.type%in%"archetype") tmp <- matrix(row.exp, nrow = nrow(X), ncol = G)
    if(prediction.type%in%"species") tmp <- matrix(row.exp, nrow = nrow(X), ncol = S)
    bPreds$fit <- tmp
    tmp.grp <- sweep(bootPreds, 1, row.exp, "-")
    tmp.grp <- tmp.grp^2
    tmp.grp <- sqrt(rowSums(tmp.grp)/(nboot - 1))
    if(prediction.type%in%"archetype") tmp.grp <- matrix(tmp.grp, nrow = nrow(X), ncol = G)
    if(prediction.type%in%"species") tmp.grp <- matrix(tmp.grp, nrow = nrow(X), ncol = S)
    bPreds$ses <- tmp.grp
    if(prediction.type%in%"archetype") colnames(bPreds$fit) <- colnames(bPreds$ses) <- object$names$SAMs
    if(prediction.type%in%"species") colnames(bPreds$fit) <- colnames(bPreds$ses) <- object$names$spp
    tmp.fun <- function(x) return(quantile(bootPreds[x, ],
                                           probs = c(0, alpha) + (1 - alpha)/2,
                                           na.rm = TRUE))
    tmp1 <- parallel::mclapply(seq_len(nrow(bootPreds)), tmp.fun,
                               mc.cores = mc.cores)
    tmp1 <- do.call("rbind", tmp1)
    if(prediction.type%in%"archetype"){
      tmp1 <- array(tmp1, c(nrow(X), G, 2), dimnames = list(NULL,
                                                            NULL, NULL))
      bPreds$cis <- tmp1[, 1:G, ]
      dimnames(bPreds$cis) <- list(NULL, object$names$SAMs, c("lower", "upper"))
    }
    if(prediction.type%in%"species"){
      tmp1 <- array(tmp1, c(nrow(X), S, 2), dimnames = list(NULL,
                                                            NULL, NULL))
      bPreds$cis <- tmp1[, 1:S, ]
      dimnames(bPreds$cis) <- list(NULL, object$names$spp, c("lower", "upper"))
    }

    # dimnames(bPreds$cis) <- list(NULL, nam, c("lower", "upper"))
    ret <- list(ptPreds = ptPreds, bootPreds = bPreds$fit,
                bootSEs = bPreds$ses, bootCIs = bPreds$cis)

    tmp.fun.grp <- function(x) return(quantile(bootPreds[x, ],
                                               probs = c(0, alpha) + (1 - alpha)/2, na.rm = TRUE))
    tmp1 <- parallel::mclapply(seq_len(nrow(bootPreds)), tmp.fun.grp,
                               mc.cores = mc.cores)
    tmp1 <- do.call("rbind", tmp1)
    if(prediction.type%in%"archetype"){
      tmp1 <- array(tmp1, c(nrow(X), G, 2), dimnames = list(NULL,
                                                            NULL, NULL))
      bPreds$fit_cis <- tmp1[, 1:G, ]
      dimnames(bPreds$fit_cis) <- list(NULL, object$names$SAMs, c("lower", "upper"))
    }
    if(prediction.type%in%"species"){
      tmp1 <- array(tmp1, c(nrow(X), S, 2), dimnames = list(NULL,
                                                            NULL, NULL))
      bPreds$fit_cis <- tmp1[, 1:S, ]
      dimnames(bPreds$fit_cis) <- list(NULL, object$names$spp, c("lower", "upper"))
    }
    ret <- list(ptPreds = ptPreds, bootPreds = bPreds$fit,
                bootSEs = bPreds$ses, bootCIs = bPreds$cis)
  }
  else ret <- ptPreds
  gc()
  return(ret)
}

#'@rdname print.species_mix
#'@name print.species_mix
#'@title Print a species_mix model object.
#'@param x A model object.
#'@param \\dots Ignored
#'@export
#'@examples
#'
#'#Print information about a species_mix model
#'\donttest{
#' library(ecomix)
#' set.seed(42)
#' sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),
#' collapse = ','),")~x1+x2"))
#' sp_form <- ~ 1
#' beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,1.9),3,2,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),
#' x2=stats::rnorm(100,0,2.5))
#' dat[,-1] <- scale(dat[,-1])
#' simulated_data <- species_mix.simulate(archetype_formula = sam_form,species_formula = sp_form,
#' data = dat,beta=beta,family="bernoulli")
#' fm1 <- species_mix(archetype_formula = sam_form,species_formula = sp_form,
#' data = simulated_data, family = 'bernoulli',  nArchetypes=3)
#'print(fm1)}

"print.species_mix" <-  function (x,...){
  cat(x$titbits$family, "species_mix model\n")
  cat("\nPi\n")
  print(round(x$pi,3))
  cat("\nCoefficients\n")
  print(x$coef[-3])

}

#'@rdname species_mix.multifit
#'@name print.species_mix.multifit
#'@param x A species mix multifit object
#'@param \\dots Ignored
#'@export
#'@examples
#'
#'#Print information about a species_mix model
#'\donttest{
#'print(fmods)
#'}

"print.species_mix.multifit" <-  function (x,...){

  if(x$groupselection){

  cat("A multiple fit",x$multiple_fits[[1]][[1]]$family,"species_mix model object\n\n")
  cat("You fitted models with ",paste0(x$nArchetypes,collapse = ", "),"Archetypes.\n\nA total of",x$nstart,"random starts were fitted per Archetype.\n\n")

  bics <- matrix(NA,x$nstart,max(seq_along(x$nArchetypes)))
  for ( ii in seq_len(x$nstart)){
    for ( jj in seq_along(x$nArchetypes)){
      bics[ii,jj] <- BIC(x$multiple_fits[[jj]][[ii]])
    }
  }

  best_mod_idx <- which(bics == min(bics,na.rm=TRUE), arr.ind = TRUE)
  best_mod <- x$multiple_fits[[best_mod_idx[2]]][[best_mod_idx[1]]]

  cat("The best model based on BIC has",best_mod$G,"archetypes\n")

  print(best_mod)

  cat("You can access more elements of this model using x$multiple_fits[[",best_mod_idx[2],"]][[",best_mod_idx[1],"]]\n")

  } else {

  }

}

#' @title Estimate residuals for a species_mix object
#' @rdname residuals.species_mix
#' @name residuals.species_mix
#' @param object A returned species_mix model object.
#' @param \dots additional calls for residual function
#' @param type The type of residuals to estimate. Default is "RQR" (Random Quantile Residuals).
#' But you can also simulate many Random Quantile Residuals using "SimRQR".
#' @param quiet Print out residual issues.
#' @export
#' @description  The randomised quantile residuals ("RQR", from Dunn and Smyth, 1996) are defined by
#' their marginal distribution function (marginality is over other species observations within that site;
#' see Woolley et al, in prep).

"residuals.species_mix" <- function( object, ..., type="RQR", quiet=FALSE) {
  if( ! type %in% c("RQR","deviance","pearson"))
    stop( "Unknown type of residual requested.\n Only deviance and RQR (for randomised quantile residuals) are implemented\n")
  if(type=="pearson"){
    stop("pearson residuals not implemented yet.")
  }
  if( type=="RQR"){
    resids <- matrix( NA, nrow=object$n, ncol=object$S)
    switch( object$family,
            bernoulli = { fn <- function(y,mu,logtheta) pbinom( q=y, size=1, prob=mu, lower.tail=TRUE)},
            poisson = { fn <- function(y,mu,logtheta) ppois( q=y, lambda=mu, lower.tail=TRUE)},
            ippm = { fn <- function(y,mu,logtheta) ppois( q=y, lambda=mu, lower.tail=TRUE)},
            negative.binomial = { fn <- function(y,mu,logtheta) pnbinom( q=y, mu=mu, size=1/exp( logtheta), lower.tail=TRUE)},
            gaussian = { fn <- function(y,mu,logtheta) pnorm( q=y, mean=mu, sd=exp( logtheta), lower.tail=TRUE)})


    for( ss in 1:object$S){
      if( object$family %in% c("bernoulli","poisson","ippm","negative.binomial")){
        tmpLower <- fn( object$titbits$Y[,ss]-1, object$mus[,ss,], object$coef$theta[ss])
        tmpUpper <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$theta[ss])
        tmpLower <- rowSums( tmpLower * object$pi)
        tmpLower <- ifelse( tmpLower<0, 0, tmpLower) #get rid of numerical errors for really small negative values
        tmpLower <- ifelse( tmpLower>1, 1, tmpLower) #get rid of numerical errors for 1+epsilon.
        tmpUpper <- rowSums( tmpUpper * object$pi)
        tmpUpper <- ifelse( tmpUpper<0, 0, tmpUpper) #get rid of numerical errors for really small negative values
        tmpUpper <- ifelse( tmpUpper>1, 1, tmpUpper) #get rid of numerical errors for 1+epsilon.
        resids[,ss] <- runif( object$n, min=tmpLower, max=tmpUpper)
        resids[,ss] <- qnorm( resids[,ss])
      }
      if( object$family == "gaussian"){
        tmp <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$theta[ss])
        tmp <- rowSums( tmp * object$pi)
        resids[,ss] <- qnorm( tmp)
      }
    }
    if(!quiet & sum( resids==Inf | resids==-Inf | is.na(resids))>0)
      message( "Some residuals, well ",sum( resids==Inf | resids==-Inf| is.na(resids)), " of ",object$n*object$S," observations are very large (infinite actually).\nThese observations lie right on the edge of the realistic range of the model for the data (maybe even over the edge).")

  }
  return( resids)
}

#' @rdname summary.species_mix
#' @name summary.species_mix
#' @title Print a summary of the species_mix model
#' @description Print a summary of the species_mix model, this requires the variance-covariance matrix to be appended to the species_mix model.
#' @param object A species_mix model object
#' @param \\dots Ignored
#' @export
"summary.species_mix" <-function (object, ...){
  if (is.null(object$vcov)) {
    object$vcov <- matrix(NA, nrow = length(unlist(object$coef)),
                          ncol = length(unlist(object$coef)))
    stop("No variance matrix has been supplied")

  }
  res <- cbind(unlist(object$coefs),  sqrt(diag(object$vcov)))
  res <- cbind(res, res[, 1]/res[, 2])
  res <- cbind(res, 2 * (1 - pnorm(abs(res[, 3]))))
  colnames(res) <- c("Estimate", "SE", "z-score", "p")
  return(res)
}

#' @rdname terms.species_mix
#' @name terms.species_mix
#' @title Return terms from a fitted species_mix model
#' @description This function returns a list of terms objects, one for each formula in the species mix model.
#' \describe{
#'  \item{xterms}{Are the terms for the archetypes formula.}
#'  \item{wterms}{Are the terms for the species formula.}
#'  \item{uterms}{Are the terms for the all/bias formula.}
#' }
#' @param x A species_mix model object
#' @param \\dots Ignored
#' @export
"terms.species_mix" <-function (x, ...){
  tt <- x$terms
  return(tt)
}



#'@rdname vcov.species_mix
#'@name vcov.species_mix
#'@title Estimate the Variance-covariance matrix for a species_mix object.
#'@description Calculates variance-covariance matrix from a species_mix object
#'@param object an object obtained from fitting a RCP (for region of common profile) mixture model. Such as that generated from a call to species_mix(qv).
#'@param object2 an object of class \code{species_mix} containing bootstrap samples of the parameter estimates (see species_mix_boot(qv)). If NULL (default) the bootstrapping is performed from within the vcov function. If not null, then the vcov estimate is obtained from these bootstrap samples.
#'@param method the method to calculate the variance-covariance matrix. Options are:'FiniteDifference' (default), \code{BayesBoot}, \code{SimpleBoot}, and \code{EmpiricalInfo}. The two bootstrap methods (\code{BayesBoot} and \code{SimpleBoot}, see species_mix_boot(qv)) should be more general and may possibly be more robust. The \code{EmpiricalInfo} method implements an empirical estimate of the Fisher information matrix, I can not recommend it however. It seems to behave poorly, even in well behaved simulations. It is computationally thrifty though.
#'@param nboot the number of bootstrap samples to take for the bootstrap estimation. Argument is ignored if !method \%in\% c(\code{FiniteDifference},'EmpiricalInfo').
#'@param mc.cores the number of cores to distribute the calculations on. Default is 4. Set to 1 if the computer is running Windows (as it cannot handle forking -- see mclapply(qv)). Ignored if method=='EmpiricalInfo'.
#'@param \\dots Ignored
#'@details If method is \code{FiniteDifference}, then the estimates variance matrix is based on a finite difference approximation to the observed information matrix.
#'If method is either "BayesBoot" or "SimpleBoot", then the estimated variance matrix is calculated from bootstrap samples of the parameter estimates. See Foster et al (in prep) for details of how the bootstrapping is actually done, and species_mix_boot(qv) for its implementation.
#'@return A square matrix of size equal to the number of parameters. It contains the variance matrix of the parameter estimates.
#'@export

#'@export
#'
#'@examples
#'
#'# Estimate the variance-covariance matrix.
#'# This will provide estimates of uncertainty for model parameters.
#'\donttest{
#' library(ecomix)
#' set.seed(42)
#' sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),
#' collapse = ','),")~x1+x2"))
#' sp_form <- ~ 1
#' beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,1.9),3,2,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),
#' x2=stats::rnorm(100,0,2.5))
#' dat[,-1] <- scale(dat[,-1])
#' simulated_data <- species_mix.simulate(archetype_formula = sam_form,species_formula = sp_form,
#' data = dat,beta=beta,family="bernoulli")
#' fm1 <- species_mix(archetype_formula = sam_form,species_formula = sp_form,
#' data = simulated_data, family = 'bernoulli',  nArchetypes=3)
#' vcov(fm1)}
"vcov.species_mix" <- function (object, object2=NULL, method = "BayesBoot",
                                nboot = 10, mc.cores = 1, ...){
  if( method %in% c("simple","Richardson"))
    method <- "FiniteDifference"
  if (!method %in% c("FiniteDifference", "BayesBoot", "SimpleBoot")) {
    error("Unknown method to calculate variance matrix, viable options are:
            'FiniteDifference' (numerical), 'BayesBoot' (bayesian bootstrap)
            and 'SimpleBoot' (case-resample bootstrap)'.")
    return(NULL)
  }
  X <- object$titbits$X
  W <- object$titbits$W
  U <- object$titbits$U
  offset <- object$titbits$offset
  spp_weights <- object$titbits$spp_weights
  site_spp_weights <- object$titbits$site_spp_weights
  y <- object$titbits$Y
  y_is_na <- object$titbits$y_is_na
  size <- object$titbits$size
  family <- object$titbits$family
  disty_cases <- c("bernoulli","poisson","ippm","negative.binomial","tweedie","gaussian","binomial")
  disty <- get_family_sam(disty_cases, family)
  S <- object$S
  G <- object$G
  n <- object$n
  npx <- object$npx
  npw <- object$npw
  npu <- object$npu
  control <- object$titbits$control

  # values for optimisation.
  inits <- object$coefs
  start_vals <- setup_inits_sam(inits, S, G, X, W, U, disty, return_list = TRUE)

  # parameters to optimise
  alpha <- as.numeric(start_vals$alpha)
  beta <- as.numeric(start_vals$beta)
  eta <- as.numeric(start_vals$eta)
  gamma <- as.numeric(start_vals$gamma)
  delta <- as.numeric(start_vals$delta)
  theta <- as.numeric(start_vals$theta)

  #scores
  getscores <- 1
  alpha.score <- as.numeric(rep(NA, length(alpha)))
  beta.score <- as.numeric(rep(NA, length(beta)))
  eta.score <- as.numeric(rep(NA, length(eta)))
  if( npw > 0){
    control$optiPart <- as.integer(1)
    gamma.score <- as.numeric(rep(NA, length(gamma)))
  } else {
    control$optiPart <- as.integer(0)
    gamma.score <- rep(-999999,S)
  }
  if(!is.null(U)) {
    npu <- as.integer(ncol(U))
    control$optiAll <- as.integer(1)
    delta.score <- as.numeric(matrix(NA, ncol=ncol(U)))
  } else {
    delta.score <- -99999
    control$optiAll <- as.integer(0)
    Ucpp <- matrix(1,nrow = n,ncol=1)
    npu <- as.integer(1) # a dummy variable to stop c++ issues.
  }
  if(disty%in%c(4,6)){
    control$optiDisp <- as.integer(1)
    theta.score <- as.numeric(rep(NA, length(theta)))
  }else{
    control$optiDisp <- as.integer(0)
    theta.score <- rep(-999999,S)
  }
  scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,gamma.score,delta.score,theta.score))))
  conv <- FALSE

  #model quantities
  pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, S))
  loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

  #remove finite for now, as we only need it for ippm
  if (method %in% c("FiniteDifference")) {
    grad_fun <- function(x) {
      x <- setup_inits_sam(x, S, G, X, W, U, disty, return_list = FALSE)
      #x is a vector of first order derivates to optimise using numDeriv in order to find second order derivates.
      start <- 0
      alpha <- x[start + seq_len(S)]
      start <- start + S
      beta <- x[start + seq_len((G*npx))]
      start <- start + (G*npx)
      eta <- x[start + seq_len(G - 1)]
      start <- start + (G-1)
      if(npw>0) {
        gamma <- x[start + seq_len((S*npw))]
        start <- start + (S*npw)
      } else {
        gamma <- rep(-999999,S)
        # start <- start + 1
      }
      if(npu>0) {
        delta <- x[start + seq_len(npu)]
        start <- start + npu
      } else {
        delta <- -999999
        # start <- start + 1
      }
      if(disty%in%c(4,6)){
        theta <- x[start + seq_len(S)]
      } else {
        theta <- rep(-999999,S)
      }
      #c++ call to optimise the model (needs pretty good starting values)
      tmp <- .Call("species_mix_cpp",
                   as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
                   as.numeric(as.matrix(W)), as.numeric(as.matrix(Ucpp)),
                   as.numeric(offset), as.numeric(spp_weights),
                   as.numeric(as.matrix(site_spp_weights)),
                   as.integer(as.matrix(!y_is_na)),
                   as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
                   as.integer(npw), as.integer(npu), as.integer(n),
                   as.integer(disty),as.integer(control$optiDisp),
                   as.integer(control$optiPart),as.integer(control$optiAll),
                   # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
                   as.double(alpha), as.double(beta), as.double(eta),
                   as.double(gamma), as.double(delta), as.double(theta),
                   as.double(powers),
                   # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
                   as.numeric(control$penalty.alpha),as.numeric(control$penalty.beta),
                   as.numeric(control$penalty.pi),as.numeric(control$penalty.gamma),
                   as.numeric(control$penalty.delta),
                   as.numeric(control$penalty.theta[1]),
                   as.numeric(control$penalty.theta[2]),
                   # SEXP &RalphaPen, SEXP &RbetaPen, SEXP &RpiPen,  SEXP &RgammaPen,
                   # SEXP &RdeltaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen,
                   alpha.score, beta.score, eta.score, gamma.score, delta.score,
                   theta.score, as.integer(control$getscores.cpp), as.numeric(scores),
                   # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
                   pis_out, mus, loglikeS, loglikeSG,
                   # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
                   as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
                   as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
                   as.integer(control$printparams.cpp),
                   # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
                   as.integer(0), as.integer(0), as.integer(1),
                   # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
                   PACKAGE = "ecomix")

      tmp1 <- c(alpha.score, beta.score, eta.score)
      names(tmp1) <- c(names(object$alpha),names(object$beta),names(object$eta))
      if( npw > 0){#class( object$titbits$species_formula) == "formula")
        tmp1 <- c(tmp1, gamma.score)
        names(tmp1) <- c(names(object$gamma),names(tmp1))
      }
      if(!is.null( object$titbits$all_formula)){
        tmp1 <- c(tmp1, delta.score)
        names(tmp1) <- c(names(object$delta),names(tmp1))
      }
      if(disty%in%c(4,6)){
        tmp1 <- c( tmp1, theta.score)
        names(tmp1) <- c(names(object$theta),names(tmp1))
      }
      gc()
      return(tmp1)
    }

    x.in <- unlist(object$coefs)
    hess <- numDeriv::jacobian(grad_fun, x = x.in, method="simple")

    hess.names <- c(names(object$alpha),names(object$beta),names(object$eta))
    if( npw > 0){#class( object$titbits$species_formula) == "formula")
      hess.names <- c(names(object$gamma),hess.names)
    }
    if(!is.null( object$titbits$all_formula)){
      hess.names <- c(names(object$delta),hess.names)
    }
    if(disty%in%c(4,6)){
      hess.names <- c(names(object$theta),hess.names)
    }

    colnames(hess)<-rownames(hess)<-hess.names

    vcov.mat <- try(-solve(hess))
    if( inherits( vcov.mat, 'try-error')){
      attr(vcov.mat, "hess") <- hess
      warning( "Hessian appears to be singular and its inverse (the vcov matrix) cannot be calculated\nThe Hessian is returned as an attribute of the result (for diagnostics).\nMy deepest sympathies.  You could try changing the specification of the model, increasing the penalties, or getting more data.")
    } else {
      vcov.mat <- ( vcov.mat + t(vcov.mat)) / 2 #to ensure symmetry
    }
  }
  if( method %in% c( "BayesBoot","SimpleBoot")){
    object$titbits$control$optimise <- TRUE #just in case it was turned off (see species_mix.multfit)
    if( is.null( object2))
      coefMat <- species_mix.bootstrap(object, nboot=nboot, type=method, mc.cores=mc.cores, quiet=TRUE)
    else
      coefMat <- object2
    vcov.mat <- cov( coefMat)
  }
  return(vcov.mat)
}
