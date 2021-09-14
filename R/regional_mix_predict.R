#' @rdname predict.regional_mix
#' @name predict.regional_mix
#' @title Predicts RCP probabilities at a series of sites. Confidence intervals are available too.
#' @param object an object obtained from fitting a RCP mixture model. Such as that generated from a call to regional_mix(qv).
#' @param object2 a regional_mix object obtained from bootstrapping the regional_mix object. Such as that generated from a call to regional_mix.bootstrap(qv). If not supplied, then predict.regional_mix will do parametric bootstrapping (otherwise non-parametric bootstrap).
#' @param newdata a data.frame (or something that can be coerced) containing the values of the covariates where predictions are to be made. If NULL (the default) then predictions are made at the locations of the original data.
#' @param nboot the number of parametric bootstrap samples to take for the bootstrap predictions, standard errors and confidence intervals. The default is 0, that is no bootstrapping is to be done and point predictions only are given. If object2 is not NULL, then the number of bootstrap samples is taken from that object (this argument is then ignored).
#' @param alpha a numeric within [0,1] (well [0.5,1] really) indicating the specified confidence for the confidence interval. Argument is redundant if nboot == 0.
#' @param mc.cores the number of cores to spread the computations over. Ignored if running on a Windows machine.
#' @param \\dots additional predict calls	ignorned
#' @details This function implements two separate, and quite different, bootstrapping routines. The first, attributable to Foster et al (2013), which implements a parametric bootstrap, whereby parameters are drawn from their sampling distribution (defined by the ML estimates and their asymptotic vcov matrix). Yes, the vcov function needs to be run first and stored in the the regional_mix object as $vcov. Typically, the vcov matrix is obtained using numerical derivatives, which can be slow to calculate and somewhat unstable/erratic. This was the original suggestion and has been superceeded by the non-parametric bootstrap routine. This is described in Foster et al (in prep) and bootstraps the sampling site data repeatedly, and for each bootstrap sample the model is re-estimated. Variation in the bootstrap samples is carried forward to the prediction step to guage the uncertainty.
#' The parametric bootsrap implementation of this function can take a while to run ??? it is a bootstrap function. nboot samples of the parameters are taken and then used to predict at each set of covariates defined in newdata. Quantiles of the resulting sets of bootstrap predictions are then taken. It is the last step that really takes a while. The non-parametric version of this function should not take as long as the grunt work of bootstrapping is carried out in the regional_mix.bootstrap(qv) function.
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

    # terms
    tt <- terms(object)

    form.X <- as.formula(object$titbit$rcp_formula)

    rcp.tm <- tt[[1]]
    dat.levels <- lapply(newdata,levels)
    mfx <- model.frame(rcp.tm, newdata, xlev = dat.levels)
    if (!is.null(cl <- attr(rcp.tm, "dataClasses")))
      .checkMFClasses(cl, mfx)
    dat.fac <- vapply(newdata, is.factor, logical(1L))
    contrasts.list <- lapply(dat.fac,function(x)ifelse(x==TRUE,"contr.treatment",NA))
    contrasts.list <- Filter(Negate(anyNA),contrasts.list)
    X <- model.matrix(rcp.tm, mfx, contrasts.arg = contrasts.list)

    if (class(object$titbits$species_formula) == "formula") {
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
    p.w <- ncol(W)
    # form.X <- as.formula(object$titbit$rcp_formula)
    # if (length(form.X) == 3)
      # form.X[[2]] <- NULL
    # X <- model.matrix(form.X, stats::model.frame(form.X, data = as.data.frame(newdata)))
    # if (class(object$titbits$species_formula) == "formula") {
      # W <- model.matrix(object$titbits$species_formula, stats::model.frame(object$titbits$species_formula,
                                                                           # data = as.data.frame(newdata)))
      # p.w <- ncol(W)
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
      message("Using supplied regional_mix.bootstrap object (non-parametric bootstrap)")
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
