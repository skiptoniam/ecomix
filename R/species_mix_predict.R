#'@rdname predict.species_mix
#'@name predict.species_mix
#'@title Predict a species_mix model.
#'@param object is a matrix model returned from the species_mix model.
#'@param boot.object is a species mix bootstrap object.
#'@param newdata a matrix of new observations for prediction.
#'@param offset an offset for prediction
#'@param nboot Number of bootstraps (or simulations if using IPPM) to run if no boot.object is provided.
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

"predict.species_mix" <- function(object, boot.object = NULL, newdata = NULL,
                                  offset = NULL, nboot = 0, alpha = 0.95,
                                  mc.cores = 1, type = 'response',
                                  prediction.type='archetype',
                                  na.action = "na.pass", ...){
  if (is.null(newdata)) {
    X <- object$titbits$X
    W <- object$titbits$W
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

  spp_wts <- object$titbits$spp_weights
  site_spp_wts <- object$titbits$site_spp_weights

  # disty_cases <- c("bernoulli","poisson","ippm","negative.binomial","tweedie","gaussian","binomial")
  # disty <- get_family_sam(disty_cases, object$titbits$family)
  disty <- object$disty
  tau <- object$tau

  # training data, needed to recompute tau (species-archetype posterior membership)
  # per bootstrap replicate -- tau is derived from the fitted parameters and the
  # training data, not from the (possibly different) prediction covariates above.
  Ytrain <- object$titbits$Y
  Xtrain <- object$titbits$X
  Wtrain <- object$titbits$W
  offsetTrain <- object$titbits$offset
  if (is.null(boot.object)) {
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
    allCoBoot <- as.matrix(boot.object)
    nboot <- nrow(boot.object)
  }
  if (is.null(allCoBoot))
    return(NULL)

  alphaBoot <- allCoBoot[, seq_len(S), drop=FALSE]
  betaBoot <- allCoBoot[, S + seq_len((G*npx)), drop=FALSE]
  etaBoot <- allCoBoot[, S + (G*npx) + seq_len(G-1), drop=FALSE]
  gammaBoot <- allCoBoot[, S + (G-1) + (G*npx) + seq_len((S*npw)), drop=FALSE]
  if (disty%in%c(3,4,5)) {
    thetaBoot <- allCoBoot[, S + (G-1) + (G*npx) + (S*npw) + seq_len(S), drop=FALSE]
  } else {
    thetaBoot <- NULL
  }

  outcomes <- matrix(NA, nrow = nrow(X), ncol = S)
  myContr <- object$titbits$control
  nam <- paste("G", 1:G, sep = "_")

  recompute_tau_sam <- function(alpha_ii, beta_ii, eta_ii, gamma_ii, theta_ii){
    pi_ii <- additive_logistic(eta_ii)
    fits_ii <- list(alpha = alpha_ii, beta = beta_ii, gamma = gamma_ii, theta = theta_ii)
    logls_ii <- get_logls_sam(y = Ytrain, X = Xtrain, W = Wtrain, G = G, S = S,
                              spp_weights = spp_wts, site_spp_weights = site_spp_wts,
                              offset = offsetTrain, disty = disty, linky = object$link,
                              size = object$titbits$size, powers = object$titbits$powers,
                              control = myContr, fits = fits_ii, get_fitted = FALSE)
    get_taus(pi_ii, logls_ii$logl_sp, G, S)
  }

  windowsOS <- Sys.info()["sysname"] == "Windows"
  if (windowsOS & mc.cores>1 & !object$titbits$control$quiet)
    message("Parallelised version of function not available for Windows machines. Reverting to single processor.")
  tauMcCores <- if(windowsOS) 1 else mc.cores

  # tau (species-archetype posterior membership) is derived from the training data and
  # is not fixed across bootstrap replicates -- recompute it once per replicate here, then
  # hand the whole (nboot x S x G) array to the C++ prediction routine below.
  if (nboot > 0) {
    tau_list <- parallel::mclapply(seq_len(nboot), function(ii){
      theta_ii <- if(!is.null(thetaBoot)) thetaBoot[ii,] else rep(-999999,S)
      recompute_tau_sam(alpha_ii = alphaBoot[ii,],
                        beta_ii = matrix(betaBoot[ii,],G,npx),
                        eta_ii = etaBoot[ii,],
                        gamma_ii = matrix(gammaBoot[ii,],S,npw),
                        theta_ii = theta_ii)
    }, mc.cores = tauMcCores)
    tauBootArr <- array(0, dim = c(nboot,S,G))
    for (ii in seq_len(nboot)) tauBootArr[ii,,] <- tau_list[[ii]]
  }

  predtype_int <- switch(prediction.type, archetype = 0L, species = 1L,
                         stop("prediction.type not known"))
  type_int <- switch(type, response = 0L, link = 1L, stop("type not known"))
  linky_int <- if(object$link=="cloglog") 1L else 0L
  Wcpp <- if(npw>0) W[,-1,drop=FALSE] else matrix(0,nrow(X),1)
  gammaIn <- if(npw>0) object$coefs$gamma else matrix(0,S,1)

  if (nboot==0) {
    bootalpha_arg <- bootbeta_arg <- bootgamma_arg <- boottau_arg <- as.numeric(0)
  } else {
    bootalpha_arg <- as.numeric(alphaBoot)
    bootbeta_arg <- as.numeric(betaBoot)
    bootgamma_arg <- if(npw>0) as.numeric(gammaBoot) else as.numeric(matrix(0,nboot,1))
    boottau_arg <- as.numeric(tauBootArr)
  }

  cpp_res <- .Call("sam_cpp_pred",
                   as.numeric(as.matrix(X)), as.numeric(as.matrix(Wcpp)), as.numeric(offset),
                   as.integer(G), as.integer(S), as.integer(nrow(X)), as.integer(npx), as.integer(npw),
                   as.integer(disty), as.integer(linky_int), as.integer(type_int), as.integer(predtype_int),
                   as.numeric(object$coefs$alpha), as.numeric(object$coefs$beta),
                   as.numeric(gammaIn), as.numeric(tau),
                   bootalpha_arg, bootbeta_arg, bootgamma_arg, boottau_arg,
                   as.integer(nboot),
                   PACKAGE = "ecomix")

  predCols <- if(prediction.type=="archetype") G else S
  predNames <- if(prediction.type=="archetype") object$names$SAMs else object$names$spp

  ptPreds <- matrix(cpp_res$preds, nrow = nrow(X), ncol = predCols)
  colnames(ptPreds) <- predNames

  if (nboot == 0) {
    gc()
    return(ptPreds)
  }

  bootPreds <- matrix(cpp_res$bootPreds, nrow = nrow(X) * predCols, ncol = nboot)
  bPreds <- list()
  row.exp <- rowMeans(bootPreds)
  bPreds$fit <- matrix(row.exp, nrow = nrow(X), ncol = predCols)
  tmp.grp <- sweep(bootPreds, 1, row.exp, "-")
  tmp.grp <- tmp.grp^2
  tmp.grp <- sqrt(rowSums(tmp.grp)/(nboot - 1))
  bPreds$ses <- matrix(tmp.grp, nrow = nrow(X), ncol = predCols)
  colnames(bPreds$fit) <- colnames(bPreds$ses) <- predNames

  tmp.fun <- function(x) return(quantile(bootPreds[x, ],
                                         probs = c(0, alpha) + (1 - alpha)/2,
                                         na.rm = TRUE))
  tmp1 <- parallel::mclapply(seq_len(nrow(bootPreds)), tmp.fun,
                             mc.cores = tauMcCores)
  tmp1 <- do.call("rbind", tmp1)
  tmp1 <- array(tmp1, c(nrow(X), predCols, 2), dimnames = list(NULL, NULL, NULL))
  bPreds$cis <- tmp1[, 1:predCols, ]
  dimnames(bPreds$cis) <- list(NULL, predNames, c("lower", "upper"))

  ret <- list(ptPreds = ptPreds, bootPreds = bPreds$fit,
              bootSEs = bPreds$ses, bootCIs = bPreds$cis)
  gc()
  return(ret)
}
