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

  # disty_cases <- c("bernoulli","poisson","ippm","negative.binomial","tweedie","gaussian","binomial")
  # disty <- get_family_sam(disty_cases, object$titbits$family)
  disty <- object$disty
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
