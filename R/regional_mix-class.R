##### Main species mix functions to export #####

#' @title regional_mix objects
#' @rdname regional_mix
#' @name regional_mix
#' @description creates an \code{regional_mix} model.
#' @param rcp_formula an object of class "formula" (or an object that can be coerced to that class). The response variable (left hand side of the formula) needs to be either 'presence', 'occurrence', 'abundance', 'biomass' or 'quantity' this will help specify the type of data to be modelled, if the response variable is disperate to the model distribution an error will be thrown. The dependent variables (the right hind side) of this formula specifies the dependence of the region of common profile (rcp) probabilities on covariates.
#' @param species_formula an object of class "formula" (or an object that can be coerced to that class). The left hand side of this formula should be left empty (it is removed if it is not empty). The right hand side of this formula specifies the dependence of the species"'" data on covariates (typically different covariates to \code{rcp_formula} to avoid confusing confounding). An example formula is observations ~ gear_type + time_of_day, where gear_type describes the different sampling gears and time_of_day describes the time of the sample. #maybe could call this detection/bias
#' @param data a List which contains named objects 'species_data': a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned. And 'covariate_data' a data frame containing the covariate data for each site. Names of columns must match that given in \code{rcp_formula} and \code{species_formula}.
#' @param nRCP The number of mixing components (groups) to fit.
#' @param distribution The family of statistical distribution to use within the ecomix models. a  choice between "bernoulli", "poisson", "negative_binomial" and "gaussian" distributions are possible and applicable to specific types of data.
#' @param offset a numeric vector of length nrow( data) that is included into the model as an offset. It is included into the conditional part of the model where conditioning is performed on the unobserved RCP type. Note that offsets cannot be included as part of the rcp_formula or species_formula arguments ??? only through this argument.
#' @param weights a numeric vector of length nrow( data) that is used as weights in the log-likelihood calculations. If NULL (default) then all weights are assumed to be identically 1.
#' @param control a list of control parameters for optimisation and calculation. See details. From \code{control} control.
#' @param inits a character string which defines the method used to initialise finite mixture model clustering. #Will have to synergise this function call across RCP and SpeciesMix. Looks like SpeciesMix uses a em.prefit to setup initialisations. regional_mix has a number of methods. This seems like a good place to setup the bivariate clustering step - cobra function.
#' @param titbits either a boolean or a vector of characters. If TRUE (default for regional_mix(qv)), then some objects used in the estimation of the model"'"s parameters are returned in a list entitled "titbits" in the model object. Some functions, for example plot.regional_mix(qv) and predict.regional_mix(qv), will require some or all of these pieces of information. If titbits=FALSE (default for regional_mix.multifit(qv)), then an empty list is returned. If a character vector, then just those objects are returned. Possible values are:"Y" for the outcome matrix, "X" for the model matrix for the RCP model, "W" for the model matrix for the species-specific model, "offset" for the offset in the model, "wts" for the model weights, "form.RCP" for the formula for the RCPs, "form.spp" for the formula for the species-specific model, "control" for the control arguments used in model fitting, "distribution" for the conditional distribution of the species data, and "power" for the power parameters used (only used in Tweedie models). Care needs to be taken when using titbits=TRUE in regional_mix.multifit(qv) calls as titbits is created for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then setting titbits=TRUE may give users problems with memory. It is more efficient, from a memory perspective, to refit the "best" model using regional_mix(qv) after identifying it with regional_mix.multifit(qv). See examples for illustration about how to do this.
#' @param power a numeric vector (length either 1 or the number of species) defining the power parameter to use in the Tweedie models. If length(power)==1, then the same power parameter is used for all species. If length(power)==No_species, then each species gets its own power parameter. Power values must be between 1 and 2, for computational reasons they should be well away from the boundary. The default is 1.6 as this has proved to be a good ball-park value for the fisheries data that the developer has previously analysed.
#' @importFrom graphics abline hist legend lines matplot par plot points polygon rect
#' @importFrom stats as.formula binomial cooks.distance cov cutree dbinom dist dnbinom dnorm dpois
#' fitted gaussian glm hclust lm logLik model.matrix model.offset model.response
#' model.weights pbinom pnbinom pnorm poisson
#' ppois predict qnorm qqnorm quantile rbinom
#' residuals rgamma rnbinom rnorm rpois runif
#' sd uniroot update update.formula
#' @export
#' @examples
#' \dontrun{
#' simulated_data <- regional_mix.simulate()
#' rcp_form <- as.formula(paste0("cbind(",paste(colnames(simulated_data[,1:20]),
#' collapse = ','),")~1+x1+x2+x3"))
#' spp_form <- observations ~ 1 + w1 + w2
#' data <- make_mixture_data(species_data = simulated_data$species_data,
#'                           covariate_data = simulated_data$covariate_data[,-1])
#' fm_regional_mix <- regional_mix(rcp_form,spp_form,data=data,distribution='bernoulli',n_mixtures=5)
#' }
"regional_mix" <- function (rcp_formula = NULL, species_formula = NULL, data,
                            nRCP = 3, distribution="bernoulli", offset=NULL,
                            weights=NULL, control = list(), inits="random2",
                            titbits = TRUE, power=1.6){
  #the control parameters
  control <- set.control( control)
  if( !control$quiet)
    message( "RCP modelling")
  call <- match.call()
  if( !is.null(rcp_formula))
    rcp_formula <- as.formula( rcp_formula)
  else{
    if( !control$quiet)
      message( "There is no RCP model!  Please provide a model
               (intercept at least) -- exitting now")
    return( NULL)
  }
  if( !is.null( species_formula))
    species_formula <- as.formula( species_formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.exclude"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  ##data <- as.data.frame(data)
  #get the data model frames and strip out any NAs
  dat <- clean_data_rcp( mf, rcp_formula, species_formula)
  #get the outcomes
  outcomes <- model.response(dat$mf.X)
  S <- check.outcomes1(outcomes)
  if (!S) {
    if( !control$quiet)
      message("Two species have the same name -- exitting now")
    return(NULL)
  }
  if( !control$quiet)
    message( "There are: ", nRCP, "RCPs to group the sites into")
  #get the design matrix for RCP part of model
  X <- get_X_rcp(rcp_formula, dat$mf.X)
  p.x <- ncol( X)
  #get design matrix for spp part of the model -- if there is one
  W <- get_W_rcp( species_formula, dat$mf.W)
  if( all( W != -999999))
    p.w <- ncol( W)
  else
    p.w <- 0
  #get offset (if not specified then it will be zeros)
  offy <- get_offset_rcp( mf, dat$mf.X, dat$mf.W)
  #get model wts (if not specified then it will be ones)
  wts <- get_wts_rcp( mf)
  #get distribution
  disty.cases <- c("bernoulli","poisson","negative_binomial","tweedie",
                   "gaussian")
  disty <- get_dist_rcp( disty.cases, distribution)
  #get power params for Tweedie
  power <- get_power_rcp( disty, power, S)
  #summarising data to console
  print.data.summ( data, dat, S, rcp_formula, species_formula, disty.cases,
                   disty, control$quiet)

  tmp <- regional_mix.fit( outcomes, W, X, offy, wts, disty, nRCP, power, inits,
                           control, nrow( X), S, p.x, p.w)

  tmp$distribution <- disty.cases[disty]
  #calculate the posterior probs
  if( nRCP>1)
    tmp$postProbs <- calcPostProbs( tmp$pis, tmp$logCondDens)
  else
    tmp$postProbs <- rep( 1, nrow( X))
  #Residuals --not calculating residuals here.
  #Information criteria
  tmp <- calcInfoCrit( tmp)
  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_rcp( titbits, outcomes, X, W, offy, wts,
                                  rcp_formula, species_formula, control,
                                  distribution, p.w=p.w, power)
  tmp$titbits$disty <- disty
  #the last bit of the regional_mix object puzzle
  tmp$call <- call

  gc()
  tmp <- tmp[sort( names( tmp))]

  class(tmp) <- "regional_mix"
  return(tmp)

  #documentation needs to be adjusted to fit new model.

}

"regional_mix.fit" <- function(outcomes, W, X, offy, wts, disty, nRCP, power,
                               inits, control, n, S, p.x, p.w){
    if( nRCP==1){ #if there is just one RCP type -- ie no dependence on environment
      tmp <- noRCPfit(outcomes, W, X, offy, wts, disty, nRCP,
                      power, inits, control, n, S, p.x, p.w)
      return( tmp)
    }

    #initial values
    start.vals <- get_start_vals_rcp( outcomes, W, X, offy, wts, disty,
                                      nRCP, S, power, inits, control$quiet)
    #doing the optimisation
    if( !control$quiet)
      message( "Quasi-Newton Optimisation")
    # if( disty != 4){ #not Tweedie
    optimiseDisp <- TRUE
    tmp <- notTweedieOptimise( outcomes, X, W, offy, wts, S, nRCP, p.x,
                               p.w, nrow( X), disty, start.vals, power, control)
    # }
    # else #Tweedie -- quite convoluted in comparison
    #   tmp <- TweedieOptimise( outcomes, X, W, offy, wts, S, nRCP, p.x, p.w, nrow( X), disty, start.vals, power, control)

    return( tmp)

  }

#' @rdname regional_mix
#' @name regional_mix.multifit
#' @param nstart for regional_mix.multifit only. The number of random starts to perform for re-fitting. Default is 10, which will need increasing for serious use.
#' @param mc.cores for regional_mix.multifit only. The number of cores to spread the re-fitting over.
#' @export

"regional_mix.multifit" <-  function(rcp_formula = NULL, species_formula = NULL,
                                     data, nRCP = 3, distribution="bernoulli",
                                     offset=NULL, weights=NULL,
                                     control = list(), inits = "random2",
                                     titbits = FALSE, power=1.6, nstart=10,
                                     mc.cores=1) {
    #the control parameters
    control <- set.control( control)
    if( !control$quiet)
      message( "RCP modelling")
    call <- match.call()
    if( !is.null(rcp_formula))
      rcp_formula <- as.formula( rcp_formula)
    else{
      if( !control$quiet)
        message( "There is no RCP model!  Please provide a model (intercept at least) -- exitting now")
      return( NULL)
    }
    if( !is.null( species_formula))
      species_formula <- as.formula( species_formula)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("data","offset","weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.exclude"
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ##data <- as.data.frame(data)
    #get the data model frames and strip out any NAs
    dat <- clean_data_rcp( mf, rcp_formula, species_formula)
    #get the outcomes
    outcomes <- model.response(dat$mf.X)
    S <- check.outcomes1(outcomes)
    if (!S) {
      if( !control$quiet)
        message("Two species have the same name -- exitting now")
      return(NULL)
    }
    if( !control$quiet)
      message( "There are: ", nRCP, "RCPs to group the sites into")
    #get the design matrix for RCP part of model
    X <- get_X_rcp(rcp_formula, dat$mf.X)
    p.x <- ncol( X)
    #get design matrix for spp part of the model -- if there is one
    W <- get_W_rcp( species_formula, dat$mf.W)
    if( all( W != -999999))
      p.w <- ncol( W)
    else
      p.w <- 0
    #get offset (if not specified then it will be zeros)
    offy <- get_offset_rcp( mf, dat$mf.X, dat$mf.W)
    #get model wts (if not specified then it will be ones)
    wts <- get_wts_rcp( mf)
    #get distribution
    disty.cases <- c("bernoulli","poisson",
                     "negative_binomial","tweedie","gaussian")
    disty <- get_dist_rcp( disty.cases, distribution)
    # #get power params for Tweedie
    # power <- get_power_rcp( disty, power)
    #summarising data to console
    print.data.summ( data, dat, S, rcp_formula, species_formula,
                     disty.cases, disty, control$quiet)

    tmp.fun <- function(x){
      if( !control$quiet & nstart>1)
        setTxtProgressBar(pb, x)
      tmpQuiet <- control$quiet
      control$quiet <- TRUE
      dumbOut <- capture.output( tmp <- regional_mix.fit( outcomes, W, X, offy,
                                                          wts, disty, nRCP,
                                                          power, inits, control,
                                                          nrow(X), S, p.x, p.w))
      control$quiet <- tmpQuiet
      tmp$distribution <- disty.cases[disty]
      #calculate the posterior probs
      if( nRCP>1)
        tmp$postProbs <- calcPostProbs( tmp$pis, tmp$logCondDens)
      else
        tmp$postProbs <- rep( 1, nrow( X))
      #Residuals --not calculating residuals here.  Need to call residuals.regional_mix
      #Information criteria
      tmp <- calcInfoCrit( tmp)
      #titbits object, if wanted/needed.
      tmp$titbits <- get_titbits_rcp( titbits, outcomes, X, W, offy, wts,
                                      rcp_formula, species_formula, control,
                                      distribution, p.w=p.w, power)
      tmp$titbits$disty <- disty
      #the last bit of the regional_mix object puzzle
      tmp$call <- call
      class(tmp) <- "regional_mix"
      return( tmp)
    }

    #    require( parallel)
    if( !control$quiet & nstart>1)
      pb <- txtProgressBar(min = 1, max = nstart, style = 3, char = "><(('> ")
    #Fit the model many times
    many.starts <- parallel::mclapply(1:nstart, tmp.fun, mc.cores=mc.cores)

    if( !control$quiet)
      message("")

    return(many.starts)
  }

##### S3 Class exports #####
#' @rdname regional_mix
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
#' @rdname regional_mix
#' @export
#' @importFrom stats BIC
"BIC.regional_mix" <- function (object, ...){
  p <- length(unlist(object$coefs))
  k <- log(object$n)
  star.ic <- -2 * object$logl + k * p
  return(star.ic)
}

#'@rdname regional_mix
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
#'@param ... ignored
#'@param oosSize The size of the with held partitions (out-of-sample size). Use 1 (default) for leave-one-out statistics, such as Cook's distance and leave-one-out validation.
#'@param times The number of times to perform the re-estimation (the number of leave out groups). For each 1:times a random partition of the data, of size oosSize, is taken and the model is fitted to one of the partitions. It is predicted to the other partition. The exception is when oosSize=1 and times=model$n (leave-one-out). In such cases (the default too), the observations are left out one-by-one and not randomly.
#'@param mc.cores The number of cores to spread the workload over. Default is 1. Argument is useless on Windows machines ??? see ?parallel::mclapply
#'@param quiet Should printing be suppressed? Default is no, it should not. Note that in either case, printing of the iteration trace etc is suppressed for each regional_mix fit.

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

#' @rdname regional_mix
#' @name extractAIC
#' @param fit  Fitted RCP model
#' @param scale scale parameter
#' @export

"extractAIC.regional_mix" <- function (fit, scale = 1, k = 2, ...){
  n <- fit$n
  edf <- length(unlist(stats::coef(fit)))
  if (is.null(k))
    k <- 2
  aic <- -2 * logLik(fit) + k * edf
  return(c(edf, aic))
}


#' @rdname regional_mix
#' @param x a fitted regional_mix model you wish to plot
#' @param type What type of residuals to plot? 'RQR'; random quantile residuals or 'deviance' residuals.
#' @param nsim number of simulations to run
#' @param alpha.conf The bounds of the confidence intervals.
#' @param quiet Run in quiet mode.
#' @param species Which species to plot as residuals.
#' @param fitted.species What scale to plot the residuals on?
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
      newy <- as.matrix( regional_mix.simulate( nRCP=nRCP, S=S, n=n, p.x=p.x, p.w=p.w, alpha=alpha, tau=tau, beta=beta, gamma=gamma, logDisps=disp, powers=power, X=X, W=W, offset=offy, distribution=x$distribution))
      tmp <- .Call("RCP_C", as.numeric(newy[, 1:S]), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
                   as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
                   alpha, tau, beta, gamma, disp, power,
                   as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
                   alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                   pis, mus, logCondDens, logls,
                   as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
                   as.integer( FALSE), as.integer( TRUE), as.integer( FALSE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")

      allResids[, s] <- get_residuals_rcp(logls, Y, x$distribution, x$coef, nRCP, type="deviance", powers=power, quiet=TRUE)
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

#' @rdname regional_mix
#' @name plot.regional_mix_stab
#' @param x x-axis
#' @param y y-axis
#' @param minWidth min width of cuts/binning
#' @param ncuts number of cuts/bins to make
#' @param ylimmo limit of y-axis
#' @param \\dots additional plotting calls
#' @export

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

#' @rdname regional_mix
#' @name predict.regional_mix
#' @param object an object obtained from fitting a RCP mixture model. Such as that generated from a call to regimix(qv).
#' @param object2 a regional_mix object obtained from bootstrapping the regimix object. Such as that generated from a call to regiboot(qv). If not supplied, then predict.regimix will do parametric bootstrapping (otherwise non-parametric bootstrap).
#' @param newdata a data.frame (or something that can be coerced) containing the values of the covariates where predictions are to be made. If NULL (the default) then predictions are made at the locations of the original data.
#' @param nboot the number of parametric bootstrap samples to take for the bootstrap predictions, standard errors and confidence intervals. The default is 0, that is no bootstrapping is to be done and point predictions only are given. If object2 is not NULL, then the number of bootstrap samples is taken from that object (this argument is then ignored).
#' @param alpha a numeric within [0,1] (well [0.5,1] really) indicating the specified confidence for the confidence interval. Argument is redundant if nboot == 0.
#' @param mc.cores the number of cores to spread the computations over. Ignored if running on a Windows machine.
#' @param \\dots additional predict calls	ignorned
#' @details This function implements two separate, and quite different, bootstrapping routines. The first, attributable to Foster et al (2013), which implements a parametric bootstrap, whereby parameters are drawn from their sampling distribution (defined by the ML estimates and their asymptotic vcov matrix). Yes, the vcov function needs to be run first and stored in the the regimix object as $vcov. Typically, the vcov matrix is obtained using numerical derivatives, which can be slow to calculate and somewhat unstable/erratic. This was the original suggestion and has been superceeded by the non-parametric bootstrap routine. This is described in Foster et al (in prep) and bootstraps the sampling site data repeatedly, and for each bootstrap sample the model is re-estimated. Variation in the bootstrap samples is carried forward to the prediction step to guage the uncertainty.
#' The parametric bootsrap implementation of this function can take a while to run ??? it is a bootstrap function. nboot samples of the parameters are taken and then used to predict at each set of covariates defined in newdata. Quantiles of the resulting sets of bootstrap predictions are then taken. It is the last step that really takes a while. The non-parametric version of this function should not take as long as the grunt work of bootstrapping is carried out in the regiboot(qv) function.
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

#'@rdname regional_mix
#'@export
"print.regional_mix" <- function (x, ...){
  ret <- list()
  ret$Call <- x$call
  ret$Distribution <- x$distribution
  ret$coef <- stats::coef(x)
  print( ret)
  invisible(ret)
}

#' @rdname regional_mix
#' @name residuals.regional_mix
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
    switch( object$distribution,
            bernoulli = { fn <- function(y,mu,logdisp,power) pbinom( q=y, size=1, prob=mu, lower.tail=TRUE)},
            poisson = { fn <- function(y,mu,logdisp,power) ppois( q=y, lambda=mu, lower.tail=TRUE)},
            negative_binomial = { fn <- function(y,mu,logdisp,power) pnbinom( q=y, mu=mu, size=1/exp( logdisp), lower.tail=TRUE)},
            # tweedie = { fn <- function(y,mu,logdisp,power) fishMod::pTweedie( q=y, mu=mu, phi=exp( logdisp), p=power)},#CHECK!!!
            gaussian = { fn <- function(y,mu,logdisp,power) pnorm( q=y, mean=mu, sd=exp( logdisp), lower.tail=TRUE)})

    for( ss in 1:object$S){
      if( all( object$titbits$power==-999999))  tmpPow <- NULL else tmpPow <- object$titbits$power[ss]
      if( object$distribution %in% c("bernoulli","poisson","negative_binomial")){
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
      # if( object$distribution == "tweedie"){
      #   nonzero <- object$titbits$Y[,ss]>0
      #   tmpObs <- matrix( rep( object$titbits$Y[,ss], object$nRCP), ncol=object$nRCP)
      #   tmp <- matrix( fn( as.numeric( tmpObs[nonzero,]), as.numeric( object$mus[nonzero,ss,]), object$coefs$disp[ss], object$titbits$power[ss]), ncol=object$nRCP)
      #   tmp <- rowSums( tmp * object$pis[nonzero,])
      #   resids[nonzero,ss] <- qnorm( tmp)
      #   tmp <- matrix( fn( as.numeric( tmpObs[!nonzero,]), as.numeric( object$mus[!nonzero,ss,]), object$coefs$disp[ss], object$titbits$power[ss]), ncol=object$nRCP)
      #   tmp <- rowSums( tmp * object$pis[!nonzero,])
      #   resids[!nonzero,ss] <- qnorm( runif( sum( !nonzero), min=0, max=tmp))
      # }
      if( object$distribution == "gaussian"){
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
  #     sims <- regional_mix.simulate( nRCP=object$nRCP, S=object$S, n=nsim, p.x=object$p.x, p.w=object$p.w, alpha=object$coef$alpha, tau=object$coef$tau, beta=object$coef$beta, gamma=object$coef$gamma, logDisps=object$coef$disp, powers=object$titbits$power, X=X1, W=W1, offset=object$titbits$offset,distribution=object$distribution)
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



#' @rdname regional_mix
#' @name regional_mix.simulate
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
#' @param distribution Text string. Specifies the distribution of the species data. Current options are "bernoulli" (default), "poisson", "negative_binomial", "tweedie" and "gaussian.
#' @export
#' @examples
#' \dontrun{
#' #generates synthetic data
#'set.seed( 151)
#'n <- 100
#'S <- 10
#'nRCP <- 3
#'my.dist <- "negative_binomial"
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
#' distribution=my.dist, logDisp=logDisp, offset=Offy)
#'
#' }
"regional_mix.simulate" <- function(nRCP=3, S=20, n=200, p.x=3, p.w=0,
                                    alpha=NULL, tau=NULL, beta=NULL, gamma=NULL,
                                    logDisps=NULL, powers=NULL, X=NULL, W=NULL,
                                    offset=NULL, distribution="bernoulli"){
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
  if( distribution == "negative_binomial" & (is.null( logDisps) | length( logDisps) != S)){
    message( "Random values for overdispersions")
    logDisps <- log( 1 + rgamma( n=S, shape=1, scale=0.75))
  }
  if( distribution=="gaussian" & (is.null( logDisps) | length( logDisps) != S)){
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
  if( !distribution%in%c("bernoulli","poisson","negative_binomial","tweedie","gaussian")){
    message( "Distribution not found, please choose from c('bernoulli','poisson','negative_binomial','tweedie','gaussian')")
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

  if( distribution=="bernoulli")
    mu <- inv.logit(etaMu)
  if( distribution %in% c("poisson","negative_binomial"))#,"tweedie"))
    mu <- exp( etaMu)
  if( distribution == "gaussian")
    mu <- etaMu

  fitted <- matrix( NA, nrow=n, ncol=S)
  for( ii in 1:n)
    fitted[ii,] <- mu[habis[ii], ,ii]

  if( distribution=="bernoulli")
    outcomes <- matrix(rbinom(n * S, 1, as.numeric( fitted)), nrow = n, ncol = S)
  if( distribution=="poisson")
    outcomes <- matrix(rpois(n * S, lambda=as.numeric( fitted)), nrow = n, ncol = S)
  if( distribution=="negative_binomial")
    outcomes <- matrix(rnbinom(n * S, mu=as.numeric( fitted), size=1/rep(exp( logDisps), each=n)), nrow = n, ncol = S)
  if( distribution=="gaussian")
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


#'@rdname regional_mix
#'@export
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

#' @rdname regional_mix
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

#'@rdname regional_mix
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

#'@title What is the average species membership per RCP?
#'@rdname regional_mix
#'@name regional_mix.species_membership
#'@param object A RCP model
#'@param boot_object A RCP model bootstrap object
#'@param CI The confidence intervals to report the range of
#'values form bootstrap
#'@export
#'@description Extracts the average species' membership for each RCP.

"regional_mix.species_membership" <- function(object, object2=NULL, CI=c(0.025,0.975), ...){

  if(is.null(object2)){
    if(check_if_sampling(object)) type <- "single_no_sp_results"
    else type <- "single_results"
  } else {
    type <- "bootstrap_results"
  }
  partial_mus <- switch(type,
                        single_no_sp_results = partial_mus_no_species_form(object),
                        single_results = partial_mus_with_species_form(object),
                        bootstrap_results = partial_mus_from_boostrap(object, object2, CI = CI))

  return(partial_mus)
}


#### Non S3 Class objects ####
#'@title What is the average species membership per RCP?
#'@rdname regional_mix
#'@name regional_mix.species_membership
#'@param object A RCP model
#'@param boot_object A RCP model bootstrap object
#'@param CI The confidence intervals to report the range of
#'values form bootstrap
#'@export
#'@description Extracts the average species' membership for each RCP.

"regional_mix.species_membership" <- function(object, object2=NULL, CI=c(0.025,0.975), ...){

  if(is.null(object2)){
    if(check_if_sampling(object)) type <- "single_no_sp_results"
    else type <- "single_results"
  } else {
    type <- "bootstrap_results"
  }
  partial_mus <- switch(type,
                       single_no_sp_results = partial_mus_no_species_form(object),
                       single_results = partial_mus_with_species_form(object),
                       bootstrap_results = partial_mus_from_boostrap(object, object2, CI = CI))

  return(partial_mus)
}

check_if_sampling <-function(object){object$p.w<1}

partial_mus_no_species_form <- function(object, ...){

  ## what are the species taus?
  tau <- coef(object)$tau
  tau <- rbind(tau, -colSums( tau))

  ## what was the the model offset?
  offy <- object$titbits$offset

  ## what is the linear predictor (eta)
  eta <- sweep(tau, 2, object$coefs$alpha, "+") + mean(offy)

  ## what is the link function of appropriate distribution?
  if(object$distribution=="bernoulli")
    link.fun <- stats::make.link('logit')
  if(object$distribution%in%c("poisson","negative_binomial"))
    link.fun <- stats::make.link('log')
  if(object$distribution=='guassian')
    link.fun <- stats::make.link('identity')

  ## what are the partial mus: dim[nRCPs,nSpp]
  partial_mus <- link.fun$linkinv(eta)
  dimnames(partial_mus)[[2]] <- object$names$spp
  dimnames(partial_mus)[[1]] <- object$names$RCPs

  ## return the partial mus if their is no sampling artifacts (species formula).
  return(partial_mus)
}


partial_mus_with_species_form <- function(object, ... ){

  ## what are the taus?
  tau <- coef(object)$tau
  tau <- rbind(tau, -colSums( tau))

  ## what was the the model offset?
  offy <- object$titbits$offset

  ## what is the linear predictor (eta)?
  eta <- sweep(tau, 2, coef(object)$alpha, "+") + mean(offy)

  ## what is the link function of appropriate distribution?
  if(object$distribution=="bernoulli")
    link.fun <- stats::make.link('logit')
  if(object$distribution%in%c("poisson","negative_binomial"))
    link.fun <- stats::make.link('log')
  if(object$distribution=='guassian')
    link.fun <- stats::make.link('identity')

  ## probabilities for all other levels of same sampling var
  res<- lapply(seq_along(object$names$Wvars),function(jj){
    new_eta <- sweep(eta, 2, coef(object)$gamma[,object$names$Wvars[jj]], "+");
    part_mu <- link.fun$linkinv(new_eta);
    return(part_mu)})

  names(res)<- object$names$Wvars
  return(res)
}

partial_mus_from_boostrap  <- function(object, object2, CI=c(0.025,0.975)){

  #set up coefficient extraction
  taus<-grepl("tau",dimnames(object2)[[2]])
  alphas<-grepl("alpha",dimnames(object2)[[2]])

  ## what is the link function of appropriate distribution?
  if(object$distribution=="bernoulli") link.fun <- stats::make.link('logit')
  if(object$distribution%in%c("poisson","negative_binomial")) link.fun <- stats::make.link('log')
  if(object$distribution=='negative_binomial') link.fun <- stats::make.link('log')
  if(object$distribution=='guassian') link.fun <- stats::make.link('identity')

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
      part_mu <- link.fun$linkinv(tmp_eta);
      res_all[[i]]<-as.matrix(part_mu)
    }

    overall_temp<-array(unlist(res_all), dim=c( length(object$names$RCPs),length(object$names$spp),nrow(object2)))
    overall_res<-list( mean=round(apply(overall_temp, c(1,2), mean),3),
                       sd= round(apply(overall_temp, c(1,2), sd),3),
                       lower= round(apply(overall_temp, c(1,2), function(x) quantile(x, probs=CI[1])),3),
                       upper= round(apply(overall_temp, c(1,2), function(x) quantile(x, probs=CI[2])),3))

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
        part_mu <- link.fun$linkinv(new_eta);
        return(part_mu)})

      names(res)<-object$names$Wvars
      res_all[[i]] <- res
    }

    #Compile list of summaries at the sampling factor level
    samp_res <- rep(list(list()), length(object$names$Wvars))
    names(samp_res) <- object$names$Wvars

    for(k in seq_along(object$names$Wvars)){
        samp_res[[k]]<-list(mean=round(apply(simplify2array(res_all[[k]]), c(1,2), mean),3),
                            sd=round(apply(simplify2array(res_all[[k]]), c(1,2), sd),3),
                            lower=round(apply(simplify2array(res_all[[k]]), c(1,2), function(x) quantile(x, probs=CI[1])),3),
                            upper=round(apply(simplify2array(res_all[[k]]), c(1,2), function(x) quantile(x, probs=CI[2])),3))
    }

    overall_temp<-list()
    for(i in seq_len(dim(object2)[1])){
        get_vals <- res_all[[i]]
        overall_temp[[i]]<-apply(simplify2array(get_vals), c(1,2), mean)
    }

    overall_samp <-list(mean=round(apply(simplify2array(overall_temp), c(1,2), mean),3),
                   sd= round(apply(simplify2array(overall_temp), c(1,2), sd),3),
                   lower= round(apply(simplify2array(overall_temp), c(1,2), function(x) quantile(x, probs=CI[1])),3),
                   upper= round(apply(simplify2array(overall_temp), c(1,2), function(x) quantile(x, probs=CI[2])),3))
    samp_res$overall<-overall_samp
    return(samp_res)
  }
}

"calcInfoCrit" <- function( ret){
  k <- length(unlist(ret$coefs))
  ret$BIC <- -2 * ret$logl + log(ret$n) * k
  ret$AIC <- -2 * ret$logl + 2 * k
#  entro <- ret$postProbs * log( ret$postProbs)
#  EN <- -sum(entro)
#  ret$ICL <- ret$BIC + 2 * EN

  return( ret)
}


"calcPostProbs" <- function( pis, logCondDens){
  logPostProbs <- log( pis) + logCondDens #would be better to be working with log(pis) previously but c'est la vie
  mset <- apply( logPostProbs, 1, max)
  logSums <- mset + log( rowSums( exp( logPostProbs-mset)))
  logPostProbs <- logPostProbs - logSums
  postProbs <- exp( logPostProbs)

  return( postProbs)

}


"check.outcomes1" <- function( outs){
  nam <- colnames( outs)
  if( length( nam) == length( unique( nam)))
    return( length( nam))
  else
    return( FALSE)

}


"clean_data_rcp" <- function( data, form1, form2){
    mf.X <- stats::model.frame(form1, data = data, na.action = stats::na.exclude)
    if( !is.null( form2)){
      mf.W <- stats::model.frame(form2, data = data, na.action = stats::na.exclude)
      ids <- c( rownames( mf.W), rownames( mf.X))[duplicated( c( rownames( mf.W), rownames( mf.X)))]  #those rows of data that are good for both parts of the model.
      mf.X <- mf.X[rownames( mf.X) %in% ids,, drop=FALSE]
      mf.W <- mf.W[rownames( mf.W) %in% ids,, drop=FALSE]
    }
    else{
      mf.W <- NULL
      ids <- rownames( mf.X)
    }
    res <- list(ids=ids, mf.X=mf.X, mf.W=mf.W)

    return( res)
}

"get_dist_rcp" <- function( disty.cases, dist1){
  error.msg <- paste( c( "Distribution not implemented. Options are: ", disty.cases, "-- Exitting Now"), collapse=" ")
  disty <- switch( dist1, "bernoulli" = 1,"poisson" = 2,"negative_binomial" = 3,"tweedie" = 4,"gaussian" = 5,{stop( error.msg)} )
  return( disty)
}

"get_long_names_rcp" <- function( object){
#function to get the names of columns for the vcov matrix or the regional_mix_boot matrix

  #defining the column names...  Trickier than you might expect
  coef.obj <- stats::coef( object)
  colnammy <- paste( names( coef.obj$alpha), "alpha", sep="_")
  if( "tau" %in% names( coef.obj))
    colnammy <- c( colnammy, paste( paste( rep( colnames( coef.obj$tau), each=nrow( coef.obj$tau)), paste( "tau", seq_len(nrow( coef.obj$tau)), sep="_"), sep="_")))
  if( "beta" %in% names( coef.obj))
    colnammy <- c( colnammy, paste( paste( rep( colnames( coef.obj$beta), each=nrow( coef.obj$beta)), paste( "beta", seq_len(nrow(coef.obj$beta)), sep="_"), sep="_")))
  if( "gamma" %in% names( coef.obj))
    colnammy <- c( colnammy, paste(  paste( rep( rownames( coef.obj$gamma), times=ncol( coef.obj$gamma)), "gamma", sep="_"), rep( colnames( coef.obj$gamma), each=nrow( coef.obj$gamma)), sep="_"))
  if( "logDisp" %in% names( coef.obj))
    colnammy <- c( colnammy, paste( names( coef.obj$logDisp), "logDisp", sep="_"))

  return( colnammy)

}


"get_offset_rcp" <-
function( mf, mf.X, mf.W)
{
  offy <- model.offset( mf)
  if( any( offy!=0))
    return( offy)
  offy <- rep( 0, nrow( mf.X))
  return( offy)
}


"get_power_rcp" <-
function( disty, power, S)
{
  if( disty == 4){
    if( length( power) == 1)
      power <- rep(power, S)
    if( length( power) != S)
      stop( "Power parameter(s) not properly specified, exitting now")
  }
  else
    power <- -999999
  return( power)

}


"get_residuals_rcp" <-
function( site.logls, outcomes, distribution, coef, nRCP, type="deviance", powers=NULL, quiet=FALSE, nsim=1000, X, W, offy)
{
  if( ! type %in% c("deviance","RQR"))
    stop( "Unknown type of residual requested. Only deviance and RQR (for randomised quantile residuals) are implemented\n")

  if( type=="deviance"){
    resids <- sqrt( -2*site.logls)
    if( !quiet){
      message( "The sign of the deviance residuals is unknown -- what does sign mean for multiple species? Their mean is also unknown -- what is a saturated model in a mixture model?")
      message( "This is not a problem if you are just looking for an over-all fit diagnostic using simulation envelopes (cf normal and half normal plots).")
      message( "It is a problem however, when you try to see how residuals vary with covariates etc.. but the meaning of these plots needs to be considered carefully as the residuals are for multiple species anyway.")
    }
  }
  if( type=="RQR"){
    ii <- 1
    X1 <- kronecker( rep( 1, nsim), X[ii,])
    W1 <- kronecker( rep( 1, nsim), W[ii,])
    sims <- regional_mix.simulate( nRCP=nRCP, S=length( coef$alpha), n=n.sim, p.x=ncol( X), p.w=ncol( W), alpha=coef$alpha, tau=coef$tau, beta=coef$beta, gamma=coef$gamma, logDisps=coef$disp, powers=pwers, X=X1, W=W1, offset=offy,distribution=distribution)

  }

  return( resids)

}


"get_start_vals_rcp" <-
function( outcomes, W, X, offy, wts, disty, G, S, power, inits, quiet=FALSE)
{
  if( !quiet)
    message( "Obtaining starting values...")

  alpha <- rep( -999999, S)
  tau <- matrix( -999999, nrow=G, ncol=S)
  if( length( W) != 1 & !is.null( W))
    gamma <- matrix( -999999, nrow=S, ncol=ncol( W))
  else
    gamma <- -999999
  beta <- matrix( -999999, nrow=G-1, ncol=ncol( X))
  if( disty>2) disp <- rep( -999999, S) else disp <- -999999

  if (inits[1] %in% c("random","random2","hclust","noPreClust")) {
    if( !inits[1] %in% "noPreClust"){
      tmp <- dist(outcomes, method = "manhattan")
      tmp1 <- hclust(tmp, method = "ward.D2")
      tmpGrp <- cutree(tmp1, G)
      tmpX <- model.matrix( ~-1+as.factor( tmpGrp))
    }
    else{
      tmpX <- scotts.rdirichlet( n=nrow( X), alpha=rep( 5, G))
    }
    if( length( W) != 1)
      df <- cbind( tmpX, W)
    else
      df <- tmpX

    lambda.seq <- sort( unique( c( seq( from=1/0.001, to=1, length=25), seq( from=1/0.1, to=1, length=10))), decreasing=TRUE)#1/seq( from=0.001, to=1, length=100)
    if( disty == 1)
      fam <- "binomial"
    if( disty == 2 | disty == 3)
      fam <- "poisson"
    if( disty == 5)
      fam <- "gaussian"
#    if( length( W) != 1)
#      df <- cbind( model.matrix( ~-1+as.factor(tmpGrp)), W)
#    else
#      df <- cbind( model.matrix( ~-1+as.factor(tmpGrp)))
    for( ss in seq_len(ncol(outcomes))){
      if( disty != 4){
        tmp.fm <- glmnet::glmnet(y=outcomes[,ss], x=df, family=fam, offset=offy, weights=wts,
          alpha=0, #ridge penalty
          lambda=lambda.seq, #the range of penalties, note that only one will be used
          standardize=FALSE,  #don't standardize the covariates (they are already standardised)
          intercept=FALSE)  #don't give me an intercept
        locat.s <- 1/1
        my.coefs <- glmnet::coef.glmnet( tmp.fm, s=locat.s)
        if( any( is.na( my.coefs))){  #just in case the model is so badly posed that mild penalisation doesn't work...
          my.coefs <- glmnet::coef.glmnet( tmp.fm, s=lambda.seq)
          lastID <- apply( my.coefs, 2, function(x) !any( is.na( x)))
          lastID <- tail( (seq_along( lastID))[lastID], 1)
          my.coefs <- my.coefs[,lastID]
        }
      }
      else{ #Tweedie needs an unconstrained fit.  May cause problems in some cases, especially if there is quasi-separation...
        df3 <- as.data.frame( cbind( y=outcomes[,ss], offy=offy, df))
        colnames( df3)[-(1:2)] <- c( paste( "grp", 1:G, sep=""), paste( "w",seq_len(ncol(W)), sep=""))
        tmp.fm1 <- fishMod::tglm( y~-1+.-offy+offset( offy), wts=wts, data=df3, p=power[ss], vcov=FALSE, residuals=FALSE, trace=0)
        my.coefs <- c( NA, tmp.fm1$coef)
        disp[ss] <- log( tmp.fm1$coef["phi"])
        my.coefs <- my.coefs[names( my.coefs) != "phi"]
      }
      alpha[ss] <- mean( my.coefs[1+1:G])
      tau[,ss] <- my.coefs[1+1:G] - alpha[ss]
      if( length( W) != 1)
        gamma[ss,] <- my.coefs[-(1:(G+1))]
      if( disty == 3){
        tmp <- MASS::theta.mm( outcomes[,ss], as.numeric(predict( tmp.fm, s=locat.s,
                                                                   type="response",
                                                                   newx=df, newoffset=offy)),
                               weights=wts, dfr=nrow(outcomes), eps=1e-4)
        if( tmp>2)
          tmp <- 2
        disp[ss] <- log( 1/tmp)
      }
      if( disty == 5){
        preds <- as.numeric( predict(tmp.fm, s=locat.s, type="link", newx=df, newoffset=offy))
        disp[ss] <- log( sqrt( sum((outcomes[,ss] - preds)^2)/nrow( outcomes)))  #should be something like the resid standard Deviation.
      }
    }
  }
  tau <- tau[1:(G-1),]  #get rid of redundant parmaeters

  #beta stuff in here...
  beta <- matrix(0, ncol = ncol(X), nrow = G - 1)
  #this code is a nice idea but glmnet uses a softmax link function, not an additive logistic...
#  tmp.fm <- glmnet( y=as.factor( tmpGrp), x=X[,-1], family="multinomial", alpha=0, lambda=1/seq( from=1,to=10,length=25), standardize=FALSE, intercept=TRUE)
#  beta <- t( sapply(  stats::coef( tmp.fm, s=locat.s), as.numeric))
#  beta <- beta[1:(G-1),]

  #################
  #### Important magic number
  my.sd <- mult <- 0.3

  if (inits[1] == "hclust" & !quiet)
    message( "Obtaining initial values for species' model from clustering algorithm -- no random component")
  if (inits[1] == "random") {
#    my.sd <- 0.1
    alpha <- alpha + rnorm(S, sd = my.sd)
    tau <- tau + as.numeric(matrix(rnorm((G - 1) * S, sd = my.sd), ncol = G - 1))
    beta <- beta + as.numeric(matrix(rnorm((G - 1) * ncol(X), mean = 0, sd = my.sd), ncol = ncol(X), nrow = G - 1))
    if( length( W) != 1 & !is.null( W))
      gamma <- gamma + as.numeric( matrix( rnorm( S*ncol(W), mean=0, my.sd), ncol=ncol( W), nrow=S))
    if( disty > 2)
      disp <- disp + rnorm( S, sd=my.sd)
  }
  if (inits[1] == "random2") {
    my.sd <- mult*sd( alpha); if( is.na( my.sd)) my.sd <- 0.1
    alpha <- alpha + rnorm(S, sd = my.sd)
    my.sd <- mult*sd( tau); if( is.na( my.sd)) my.sd <- 0.1
    tau <- tau + as.numeric(matrix(rnorm((G - 1) * S, sd = my.sd), ncol = G - 1))
    my.sd <- mult*sd( beta); if( is.na( my.sd) | my.sd==0) my.sd <- 0.1
    beta <- beta + as.numeric(matrix(rnorm((G - 1) * ncol(X), mean = 0, sd = my.sd), ncol = ncol(X), nrow = G - 1))
    if( length( W) != 1 & !is.null( W)){
      my.sd <- mult*sd( gamma); if( is.na( my.sd) | my.sd==0) my.sd <- 0.1
      gamma <- gamma + as.numeric( matrix( rnorm( S*ncol(W), mean=0, my.sd), ncol=ncol( W), nrow=S))
    }
    if( disty > 2){
      my.sd <- mult*sd( disp); if( is.na( my.sd) | my.sd==0) my.sd <- 0.1
      disp <- disp + as.numeric( rnorm( S, mean=0, my.sd))
#      message( "My Starting Dispersions Are: ", disp,"\n")
    }
  }
  if( any( alpha == -999999)) {
    if( !quiet)
      message("Using supplied initial values (unchecked). Responsibility is entirely the users!")
    start <- 0
    alpha <- inits[start+1:S]
    start <- start + S
    tau <- inits[start + 1:((G - 1) * S)]
    start <- start + (G-1)*S
    beta <- inits[start + 1:((G - 1) * ncol(X))]
    start <- start + (G-1)*ncol(X)
    if( length( W) != 1 & !is.null( W)){
      gamma <- inits[start+ 1:(S*ncol(W))]
      start <- start + S*ncol(W)
    }
    if( disty %in% 3:5)
      disp <- inits[start+1:S]
  }

  res <- list()
  res$alpha <- as.numeric( alpha)
  res$tau <- as.numeric( tau)
  res$beta <- as.numeric( beta)
  res$gamma <- as.numeric( gamma)
  res$disp <- as.numeric( disp)

  return( res)
}


"get_titbits_rcp" <-
function( titbits, outcomes, X, W, offset, wts, rcp_formula, species_formula, control, distribution, p.w, power)
{
  if( titbits==TRUE)
    titbits <- list( Y = outcomes, X = X, W = W, offset = offset, wts=wts, rcp_formula = rcp_formula, species_formula = species_formula, control = control, distribution = distribution, power=power)
  else{
    titbits <- list()
    if( "Y" %in% titbits)
      titbits$Y <- outcomes
    if( "X" %in% titbits)
      titbits$X <- X
    if( "W" %in% titbits)
      titbits$W <- W
    if( "offset" %in% titbits)
      titbits$offset <- offset
    if( "wts" %in% titbits)
      titbits$wts <- wts
    if( "rcp_formula" %in% titbits)
      titbits$rcp_formula <- rcp_formula
    if( "species_formula" %in% titbits)
      titbits$species_formula <- species_formula
    if( "control" %in% titbits)
      titbits$control <- control
    if( "distribution" %in% titbits)
      titbits$distribution <- distribution
    if( "power" %in% titbits)
      titbits$power <- power
  }
  if( p.w==0 & "W" %in% names( titbits))
    titbits$W <- NULL
  if( p.w!=0 & "species_formula" %in% names( titbits))
    environment( titbits$species_formula) <- environment( titbits$rcp_formula)
  return( titbits)

}


"get_W_rcp" <- function( species_formula, mf.W){
  form.W <- species_formula
  if( !is.null( species_formula)){
    if( length( form.W)>2)
      form.W[[2]] <- NULL #get rid of outcomes
    W <- model.matrix( form.W, mf.W)
    tmp.fun <- function(x){ all( x==1)}
    intercepts <- apply( W, 2, tmp.fun)
    W <- W[,!intercepts,drop=FALSE]
  }
  else
    W <- -999999
  return( W)
}


"get_wts_rcp" <-function ( mf){
  wts <- model.weights( mf)
  if( is.null( wts))
    return( rep( 1, nrow( mf)))  #all weights assumed equal
  return( wts)
}

"get_X_rcp" <-function( rcp_formula, mf.X){
  form.X <- rcp_formula
  form.X[[2]] <- NULL
  form.X <- as.formula(form.X)
  X <- model.matrix(form.X, mf.X)

  return( X)
}

"globCIFinder" <-function( x, en, alpha, nsim){
	#this now works for both upper and lower CIs
	c <- uniroot( f=globErrorFn, interval=c(0.1,5), x=x, en=en, alpha=alpha, nsim=nsim)$root
	return( en*c)

}

"globErrorFn" <-function( c1, x, en, alpha, nsim){
	if( alpha > 0.5){
		tmp <- apply( x, 2, function(x) any( x-c1*en > 0))
		return( sum( tmp) / nsim - (1-alpha))
	}
	else{
		tmp <- apply( x, 2, function(x) any( x-c1*en < 0))
		return( sum( tmp) / nsim - alpha)
	}
}

"inv.logit" <-function(x){
	eta <- exp( x)
	mu <- eta / (1+eta)
	return(mu)
}

"logLik.regional_mix" <-function (object, ...){
    return(object$logl)
}


"my.rmvnorm" <-function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), method = c("eigen", "svd", "chol")){
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*%
            t(ev$vectors)
    }
    else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)
    retval
}


"nd2" <- function(x0, f, m=NULL, D.accur=4, eps=NULL, mc.cores=getOption("mc.cores", 4L), ...) {
# A function to compute highly accurate first-order derivatives
# Stolen (mostly) from the net and adapted / modified by Scott (scott.foster@csiro.au)
# From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
##multicore approach taken Fri June 26 2015

# x0 is the point where the derivative is to be evaluated,
# f is the function that requires differentiating
# m is output dimension of f, that is f:R^n -> R^m
#D.accur is the required accuracy of the resulting derivative. Options are 2 and 4. The 2 choice does a two point finite difference approximation and the 4 choice does a four point finite difference approximation.
#eps is the finite difference step size
#mc.cores is the number of cores to spread the computations over
#... other arguments to pass to f

# Report any bugs to Scott!

#  require( parallel) #for mclapply
  D.n<-length(x0)
  if (is.null(m)) {
    D.f0<-f(x0, ...)
    m<-length(D.f0) }
  if (D.accur==2) {
    D.w<-tcrossprod(rep(1,m),c(-1/2,1/2))
    D.co<-c(-1,1) }
  else {
    D.w<-tcrossprod(rep(1,m),c(1/12,-2/3,2/3,-1/12))
    D.co<-c(-2,-1,1,2) }
  D.n.c<-length(D.co)
  if( is.null( eps)) {
    macheps<-.Machine$double.eps
    D.h<-macheps^(1/3)*abs(x0)
  }
  else
    D.h <- rep( eps, D.accur)
  D.deriv<-matrix(NA,nrow=m,ncol=D.n)
  mc.fun <- function(ii){
    D.temp.f<-matrix(0,m,D.n.c)
    for (jj in 1:D.n.c) {
      D.xd<-x0+D.h[ii]*D.co[jj]*(1:D.n==ii)
      D.temp.f[,jj]<-f(D.xd, ...) }
    ret<-rowSums(D.w*D.temp.f)/D.h[ii]
    return( ret)
  }

  tmp.fun.vals <- parallel::mclapply( 1:D.n, mc.fun, mc.cores=mc.cores)
  ret <- do.call( "rbind", tmp.fun.vals)
  return( ret)


#  for (ii in 1:D.n) {
#    D.temp.f<-matrix(0,m,D.n.c)
#    for (jj in 1:D.n.c) {
#      D.xd<-x0+D.h[ii]*D.co[jj]*(1:D.n==ii)
#      D.temp.f[,jj]<-f(D.xd, ...) }
#    D.deriv[,ii]<-rowSums(D.w*D.temp.f)/D.h[ii] }
#  return( D.deriv)
}


"noRCPfit" <- function( outcomes, W, X, offy, wts, disty, nRCP, power, inits, control, n, S, p.x, p.w){
  beta <- tau <- NULL
  if( all(W==-999999)){
    W <- matrix( 1, ncol=1, nrow=n)
    gamma <- NULL
  }
  else{
    W <- cbind( 1, W)
    gamma <- matrix( NA, nrow=S, ncol=p.w)
  }
  logls <- alpha <- rep( 0, S)
  if( disty>2) disp <- rep( NA, S) else disp <- NULL
  mus <- array( NA, dim=c( n, S, nRCP))  #container for the fitted spp model
  for( ss in seq_len(ncol(outcomes))){
    if( disty == 1){ #bernoulli
      tmp.fm <- glm( cbind( outcomes[,ss], 1-outcomes[,ss]) ~ -1+W, family=binomial(), offset=offy, weights=wts)
      logls[ss] <- sum( dbinom( outcomes[,ss], size=1, prob=tmp.fm$fitted, log=TRUE))
    }
    if( disty==2){  #poisson
      tmp.fm <- glm( outcomes[,ss] ~ -1+W, family=poisson(), offset=offy, weights=wts)
      logls[ss] <- sum( dpois( outcomes[,ss], lambda=tmp.fm$fitted, log=TRUE))
    }
    if( disty==3){  #negative_binomial
      df3 <- as.data.frame( cbind( y=outcomes[,ss], offy=offy, W))
      tmp.fm <- MASS::glm.nb( y~.-1-offy+offset(offy), data=df3, weights=wts)
      logls[ss] <- sum( dnbinom( x=outcomes[,ss], size=tmp.fm$theta, mu=tmp.fm$fitted, log = TRUE))
      disp[ss] <- log( 1/tmp.fm$theta)
    }
    # if( disty==4){  #Tweedie
    #   df3 <- as.data.frame( cbind( y=outcomes[,ss], offy=offy, W))
    #   tmp.fm <- fishMod::tglm( y~.-1-offy+offset(offy), wts=wts, data=df3, p=power[ss], vcov=FALSE, residuals=FALSE, trace=0)
    #   logls[ss] <- sum( fishMod::dTweedie( y=outcomes[,ss], mu=tmp.fm$fitted, phi=tmp.fm$coef["phi"], p=power[ss], LOG=TRUE))
    #   disp[ss] <- log( tmp.fm$coef["phi"])
    #   tmp.fm$coef <- tmp.fm$coef[names( tmp.fm$coef) != "phi"]
    # }
    if( disty==5){  #gaussian
      tmp.fm <- lm( outcomes[,ss] ~-1+W, offset=offy, weights=wts)
      disp[ss] <- log(summary(tmp.fm)$sigma)
      logls[ss] <- sum( sqrt(2*pi)+exp(disp[ss])+dnorm( outcomes[,ss], mean=tmp.fm$fitted, sd=exp(disp[ss]), log=TRUE))
    }
    alpha[ss] <- tmp.fm$coef[1]
    if( p.w>0)
      gamma[ss,] <- tmp.fm$coef[-1]
    mus[,ss,1] <- fitted( tmp.fm)

  }
  logl <- sum( logls)
  #add on the penalties
  #no penalty for pi, as log(1)=0
  #no penalty for tau as all zero
  #gamma
  if( !is.null( gamma))
    logl <- logl - sum( (gamma^2)/(2*control$penalty.gamma*control$penalty.gamma))
  #disp
  if( !is.null( disp))
    logl <- logl - sum( ((disp-control$penalty.disp[1])^2)/(2*control$penalty.disp[2]*control$penalty.disp[2]))

  ret <- list()
  ret$pis <- matrix(1, ncol = nRCP, nrow=n)
  ret$mus <- array( mus, dim=c(n,S,nRCP))
  ret$coefs <- list(alpha = alpha, tau = NULL, beta = NULL, gamma=gamma, disp=disp)
  ret$scores <- NULL
  ret$logCondDens <- NULL
  ret$conv <- NULL
  ret$S <- S; ret$nRCP <- nRCP; ret$p.x <- p.x; ret$p.w <- p.w; ret$n <- n
  ret$start.vals <- NULL
  ret$logl <- logl
  ret$logl.sites <- NULL  #for residuals

  return( ret)
}


"notTweedieOptimise" <-function( outcomes, X, W, offy, wts, S, nRCP, p.x, p.w, n, disty, start.vals, power, control){
  inits <- c(start.vals$alpha, start.vals$tau, start.vals$beta, start.vals$gamma, start.vals$disp)
  alpha <- start.vals$alpha; tau <- as.numeric( start.vals$tau); beta <- as.numeric( start.vals$beta); gamma <- as.numeric( start.vals$gamma); disp <- start.vals$disp
  #scores
  alpha.score <- as.numeric(rep(NA, S))
  tau.score <- as.numeric(matrix(NA, ncol = S, nrow = nRCP - 1))
  beta.score <- as.numeric(matrix(NA, ncol = ncol(X), nrow = nRCP - 1))
  if( p.w > 0)
    gamma.score <- as.numeric(matrix( NA, nrow=S, ncol=ncol(W)))
  else
    gamma.score <- -999999
  if( disty %in% 3:5)
    disp.score <- as.numeric( rep( NA, S))
  else
    disp.score <- -999999
  scoreContri <- -999999#as.numeric(matrix(NA, ncol = length(inits), nrow = n))
  #model quantities
  pis <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c( n, S, nRCP)))  #container for the fitted spp model
  logCondDens <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))
  logls <- as.numeric(rep(NA, n))
  conv <- as.integer(0)

  tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
          as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
          alpha, tau, beta, gamma, disp, power,
          as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
          alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
          pis, mus, logCondDens, logls,
          as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
          as.integer( control$optimise), as.integer(control$loglOnly), as.integer( control$derivOnly), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")

  ret <- list()
  ret$pis <- matrix(pis, ncol = nRCP)
  ret$mus <- array( mus, dim=c(n,S,nRCP))
  ret$coefs <- list(alpha = alpha, tau = tau, beta = beta, gamma=gamma, disp=disp)
  if( any( ret$coefs$gamma==-999999, na.rm=TRUE))
    ret$coefs$gamma <- NULL
  if( any( ret$coefs$disp==-999999, na.rm=TRUE))
    ret$coefs$disp <- NULL
  ret$names <- list( spp=colnames( outcomes), RCPs=paste( "RCP", 1:nRCP, sep=""), Xvars=colnames( X))
  if( p.w>0)
    ret$names$Wvars <- colnames( W)
  else
    ret$names$Wvars <- NA
  ret$scores <- list(alpha = alpha.score, tau = tau.score, beta = beta.score, gamma = gamma.score, disp=disp.score)
  if( any( ret$scores$gamma==-999999, na.rm=TRUE))
    ret$scores$gamma <- NULL
  if( any( ret$scores$disp==-999999, na.rm=TRUE))
    ret$scores$disp <- NULL
  ret$logCondDens <- matrix(logCondDens, ncol = nRCP)
  if( control$optimise)
    ret$conv <- conv
  else
    ret$conv <- "not optimised"
  ret$S <- S; ret$nRCP <- nRCP; ret$p.x <- p.x; ret$p.w <- p.w; ret$n <- n
  ret$start.vals <- inits
  ret$logl <- tmp
  ret$logl.sites <- logls  #for residuals

  return( ret)

}

#
# "orderFitted" <- function( fm, simDat){
# 	RCPs <- attr( simDat, "RCP")
# 	posts <- fm$postProbs
#
# 	perms <- gtools::permutations( length( unique( RCPs)), length( unique( RCPs)))
# 	classErr <- rep( NA, ncol( perms))
# 	classErrRunnerUp <- classErr
# 	for( ii in seq_len(nrow(perms))){
# 		postsTMP <- posts[,perms[ii,]]
# 		postsTMP <- apply( postsTMP, 1, which.max)
# 		my.tab <- table( RCPs, postsTMP)
# 		classErr[ii] <- sum( diag( my.tab)) / sum( my.tab)
# 	}
# 	perms <- perms[which.max( classErr),]
# 	#coefs
# 	tau <- matrix( fm$coefs$tau, nrow=fm$nRCP-1, ncol=fm$S)
# 	tau <- rbind( tau, -colSums( tau))
# 	tau <- tau[perms,]
# 	beta <- matrix( fm$coefs$beta, nrow=fm$nRCP-1, ncol=fm$p.x)
# 	beta <- rbind( beta, 0)
# 	beta <- beta[perms,]
# 	beta <- beta - rep( beta[fm$nRCP,], each=fm$nRCP)
# 	fm$coefs$tau <- as.numeric( tau[-fm$nRCP,])
#   fm$coef$beta <- as.numeric( beta[-fm$nRCP,])
# 	#scores
# 	fm$scores <- NULL
# 	#pis
# 	fm$pis <- fm$pis[,perms]
# 	#postProbs
# 	fm$postProbs <- fm$postProbs[,perms]
# 	#mus
# 	fm$mus <- fm$mus[,,perms]
# 	#vcov
# 	fm$vcov <- NULL
# 	#order
# 	fm$perm <- perms
# 	#classification error
# 	fm$classErr <- max( classErr)
# 	fm$classErrRunnerUp <- max( classErr[-(which.max( classErr))])
#
# 	return( fm)
#
# }
#

# "orderPost" <- function( new.fm=NULL, fm, RCPs=NULL, sample=NULL){
# 	G1 <- G2 <- NULL
# 	if( !is.null( new.fm))
# 		G <- G1 <- new.fm$nRCP
# 	if( !is.null( RCPs))
# 		G <- G2 <- length( unique( RCPs))
# 	if( sum( !is.null( c(G1,G2))) != 1){
# 		message( "Problem with ordering -- provide new.fm *or* RCPs, but not both!")
# 		return( NULL)
# 	}
# 	perms <- gtools::permutations( G, G)
#
# 	if( !is.null( RCPs)){
# 		fm$postProbs <- matrix( 0, nrow=nrow( fm$postProbs), ncol=ncol( fm$postProbs))
# 		for( ii in 1:fm$nRCP)
# 			fm$postProbs[,ii] <- ifelse( RCPs==ii, 1, 0)
# 	}
# 	if( !is.null( sample))
# 		fm$postProbs <- fm$postProbs[sample,]
# 	classErr <- rep( NA, ncol( perms))
# 	for( ii in seq_len(nrow(perms))){
# 		my.tab <- t(fm$postProbs) %*% new.fm$postProbs[,perms[ii,]]
# 		classErr[ii] <- sum( diag( my.tab)) / sum( my.tab)
# 	}
# 	perms <- perms[which.max( classErr),]
# 	#coefs
#   alpha <- new.fm$coefs$alpha
#   gamma <- new.fm$coefs$gamma
#   disp <- new.fm$coef$disp
# 	tau <- matrix( new.fm$coefs$tau, nrow=new.fm$nRCP-1, ncol=new.fm$S)
# 	tau <- rbind( tau, -colSums( tau))
# 	tau <- tau[perms,]
#   new.fm$coefs$tau <- as.numeric( tau[-new.fm$nRCP,])
# 	beta <- matrix( new.fm$coefs$beta, nrow=new.fm$nRCP-1, ncol=new.fm$p.x)
# 	beta <- rbind( beta, 0)
# 	beta <- beta[perms,]
# 	beta <- beta - rep( beta[new.fm$nRCP,], each=3)
#   new.fm$coefs$beta <- as.numeric( beta[-new.fm$nRCP,])
# 	#scores
# 	new.fm$scores <- NULL
# 	#pis
# 	new.fm$pis <- new.fm$pis[,perms]
# 	#postProbs
# 	new.fm$postProbs <- new.fm$postProbs[,perms]
# 	#mus
# 	new.fm$mus <- new.fm$mus[,,perms]
# 	#vcov
# 	new.fm$vcov <- NULL
# 	#order
# 	new.fm$perm <- perms
# 	#classification error
# 	new.fm$classErr <- max( classErr)
# 	new.fm$classErrRunnerUp <- max( classErr[-(which.max( classErr))])
#
# 	return( new.fm)
# }

"print.data.summ" <- function( data, dat, S, rcp_formula, species_formula, disty.cases, disty, quiet=FALSE){
  if( quiet)
    return( NULL)
  n.tot <- nrow( data)
  n <- length( dat$ids)
  message("There are ", n, " fully present observations and ", n.tot, " observations in total")
  message("There are ", S, " species")
  rcp_formula[[2]] <- NULL
  message("The model for the (latent) RCP classes is: ", Reduce( "paste", deparse( rcp_formula)))
  if( !is.null( species_formula))
    message("The model for each species is: ", Reduce( "paste", deparse( species_formula)))
  else
    message("There is NO model for each species (apart from intercept(s))")
  message("The error distribution is: ", disty.cases[disty])
}

#' @rdname regional_mix
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

#' @rdname regional_mix
#' @export

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

"scotts.rdirichlet" <-
function (n, alpha)
{
  #stolen from gtools' rdirichlet
  len <- length(alpha)
  x <- matrix( rgamma( len * n, alpha), ncol = len, byrow = TRUE)
  sm <- x %*% rep(1, len)
  return( x/as.vector(sm))
}


"set.control" <-
function(control)
{
  if (!("maxit" %in% names(control)))
    control$maxit <- 500
  if( !("quiet" %in% names( control)))
    control$quiet <- FALSE
  if (!("trace" %in% names(control)))
    control$trace <- 1
  if( control$quiet)
    control$trace <- 0  #for no tracing
  if (!("nreport" %in% names(control)))
    control$nreport <- 10
  if (!("abstol" %in% names(control)))
    control$abstol <- 1e-05
  if (!("reltol" %in% names(control)))
    control$reltol <- sqrt(.Machine$double.eps)
  if (!("optimise" %in% names( control)))
    control$optimise <- TRUE
  if (!("loglOnly" %in% names(control)))
    control$loglOnly <- TRUE
  if (!("derivOnly" %in% names( control)))
    control$derivOnly <- TRUE
  if (!("penalty" %in% names(control)))
    control$penalty <- 0.01
  else
    if (control$penalty < 0) {
      message("Supplied penalty for pis is negative, reverting to the default")
      penalty <- 0.01
    }
  if (!("penalty.tau" %in% names( control)))
    control$penalty.tau <- 10
  else
  if (control$penalty.tau <= 0) {
    message("Supplied penalty for taus is negative, reverting to the default")
    control$penalty.tau <- 10
  }
  if( !("penalty.gamma" %in% names( control)))
    control$penalty.gamma <- 10
  else
    if( control$penalty.gamma <=0){
      message("Supplied penalty for gammas is negative, reverting to the default")
      control$penalty.gamma <- 10
    }
  if( !("penalty.disp" %in% names( control)))
    control$penalty.disp <- c( 10, sqrt( 10))  #the mu and sd of a log-normal
  else
    if( control$penalty.disp[2] <= 0 | length( control$penalty.disp) != 2) {
      message("Supplied penalty parameters for the dispersions is illogical, reverting to the default")
      control$penalty.disp <- c( 10, sqrt( 10))
    }

  return( control)

}


# "TweedieOptimise" <-
# function( outcomes, X, W, offy, wts, S, nRCP, p.x, p.w, n, disty, start.vals, power, control)
# {
#   Tw.phi.func <- function( phi1, spp3){
#     disp3 <- disp
#     disp3[spp3] <- phi1
#     tmp1 <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
#       as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
#       alpha, tau, beta, gamma, disp3, power,
#       as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
#       alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
#       pis, mus, logCondDens, logls,
#       as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
#       as.integer(FALSE), as.integer(TRUE), as.integer( FALSE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")
#     return( -as.numeric( tmp1))
#   }

#  Tw.phi.func.grad <- function( phi1, spp3){
#    disp3 <- disp
#    disp3[spp3] <- phi1
#    tmp.disp.score <- rep( -99999, S)
#    tmp1 <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
#      as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
#      alpha, tau, beta, gamma, disp3, power,
#      as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
#      alpha.score, tau.score, beta.score, gamma.score, tmp.disp.score, scoreContri,
#      pis, mus, logCondDens, logls,
#      as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
#      as.integer(FALSE), as.integer(FALSE), as.integer(TRUE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")
#    return( -as.numeric( tmp.disp.score[spp3]))
#  }
#
#  inits <- c(start.vals$alpha, start.vals$tau, start.vals$beta, start.vals$gamma, start.vals$disp)
#  alpha <- start.vals$alpha; tau <- as.numeric( start.vals$tau); beta <- as.numeric( start.vals$beta); gamma <- as.numeric( start.vals$gamma); disp <- start.vals$disp
#  #scores
#  alpha.score <- as.numeric(rep(NA, S))
#  tau.score <- as.numeric(matrix(NA, ncol = S, nrow = nRCP - 1))
#  beta.score <- as.numeric(matrix(NA, ncol = ncol(X), nrow = nRCP - 1))
#  if( p.w > 0)
#    gamma.score <- as.numeric(matrix( NA, nrow=S, ncol=ncol(W)))
#  else
#    gamma.score <- -999999
#  if( disty %in% 3:5)
#    disp.score <- as.numeric( rep( NA, S))
#  else
#    disp.score <- -999999
#  scoreContri <- -999999  #as.numeric(matrix(NA, ncol = length(inits), nrow = n))
#  #model quantities
#  pis <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))  #container for the fitted RCP model
#  mus <- as.numeric(array( NA, dim=c( n, S, nRCP)))  #container for the fitted spp model
#  logCondDens <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))
#  logls <- as.numeric(rep(NA, n))
#  conv <- as.integer(0)
#
#  optimiseDisp <- FALSE
#  kount <- 1
#  tmp.new <- tmp.old <- -999999
#  if( control$optimise){
#    while( (abs( abs( tmp.new - tmp.old) / ( abs( tmp.old) + control$reltol)) > control$reltol | kount==1) & (kount < 15)){
#      kount <- kount + 1
#      tmp.old <- tmp.new
#      message( "Updating Location Parameters: ", appendLF=FALSE)
#      tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
#        as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
#        alpha, tau, beta, gamma, disp, power,
#        as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
#        alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
#        pis, mus, logCondDens, logls,
#        as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
#        as.integer(control$optimise), as.integer(TRUE), as.integer( FALSE), as.integer(optimiseDisp), as.integer( FALSE), PACKAGE = "ecomix")
#      message( "Updating Dispersion Parameters: ", appendLF=FALSE)
#      for( ii in 1:S){
#        tmp1 <- nlminb( disp[ii], Tw.phi.func, Tw.phi.func.grad, spp3=ii, control=list( trace=0))
#        disp[ii] <- tmp1$par
#        message( tmp1$objective, " ")
#      }
#      message( "")
#      tmp.new <- -tmp1$objective
#    }
#  }
#  tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
#    as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
#    alpha, tau, beta, gamma, disp, power,
#    as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
#    alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
#    pis, mus, logCondDens, logls,
#    as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
#    as.integer(FALSE), as.integer( TRUE), as.integer(TRUE), as.integer(TRUE), as.integer( FALSE), PACKAGE = "ecomix")
#
#  ret <- list()
#
#  ret$pis <- matrix(pis, ncol = nRCP)
#    ret$mus <- array( mus, dim=c(n,S,nRCP))
#  ret$coefs <- list(alpha = alpha, tau = tau, beta = beta, gamma=gamma, disp=disp)
#  if( any( ret$coefs$gamma==-999999, na.rm=TRUE))
#    ret$coefs$gamma <- NULL
#  if( any( ret$coefs$disp==-999999, na.rm=TRUE))
#    ret$coefs$disp <- NULL
#  ret$names <- list( spp=colnames( outcomes), RCPs=paste( "RCP", 1:nRCP, sep=""), Xvars=colnames( X))
#  if( p.w>0)
#    ret$names$Wvars <- colnames( W)
#  else
#    ret$names$Wvars <- NA
#  ret$scores <- list(alpha = alpha.score, tau = tau.score, beta = beta.score, gamma = gamma.score, disp=disp.score)
#  if( any( ret$scores$gamma==-999999, na.rm=TRUE))
#    ret$scores$gamma <- NULL
#  if( any( ret$scores$disp==-999999, na.rm=TRUE))
#    ret$scores$disp <- NULL
#  ret$logCondDens <- matrix(logCondDens, ncol = nRCP)
#  if( control$optimise)
#    ret$conv <- conv
#  else
#    ret$conv <- "not optimised"
#  ret$S <- S; ret$nRCP <- nRCP; ret$p.x <- p.x; ret$p.w <- p.w; ret$n <- n
#  ret$start.vals <- inits
#  ret$logl <- tmp
#  ret$logl.sites <- logls  #for residuals
#
#  return( ret)
# }

# MVB's workaround for futile CRAN 'no visible blah' check:
globalVariables( package="ecomix",
  names=c( ".Traceback"
    ,"dll.path"
    ,"libname"
    ,"pkgname"
    ,"subarch"
    ,"r_arch"
    ,"this.ext"
    ,"dynlib.ext"
    ,"dlls"
    ,"x"
    ,"tmp"
    ,"p"
    ,"object"
    ,"coefs"
    ,"k"
    ,"star.ic"
    ,"logl"
    ,"n"
    ,"ret"
    ,"logPostProbs"
    ,"pis"
    ,"logCondDens"
    ,"mset"
    ,"logSums"
    ,"postProbs"
    ,"nam"
    ,"outs"
    ,"mf.X"
    ,"form1"
    ,"data"
    ,"form2"
    ,"mf.W"
    ,"ids"
    ,"res"
    ,"alpha"
    ,"spp"
    ,"tau"
    ,"nRCP"
    ,"S"
    ,"p.x"
    ,"Xvars"
    ,"p.w"
    ,"Wvars"
    ,"disp"
    ,"logDisp"
    ,"oosSize"
    ,"model"
    ,"titbits"
    ,"quiet"
    ,"pb"
    ,"txtProgressBar"
    ,"times"
    ,"funny"
    ,"setTxtProgressBar"
    ,"OOBag"
    ,"inBag"
    ,"new.wts"
    ,"wts"
    ,"control"
    ,"tmpmodel"
    ,"Y"
    ,"W"
    ,"X"
    ,"disty"
    ,"OOSppPreds"
    ,"ss"
    ,"mus"
    ,"newPis"
    ,"r.negi"
    ,"alpha.score"
    ,"tau.score"
    ,"beta.score"
    ,"gamma.score"
    ,"disp.score"
    ,"scoreContri"
    ,"logls"
    ,"conv"
    ,"tmplogl"
    ,"penalty"
    ,"penalty.tau"
    ,"penalty.gamma"
    ,"penalty.disp"
    ,"maxit"
    ,"nreport"
    ,"abstol"
    ,"reltol"
    ,"ret.logl"
    ,"mc.cores"
    ,"parallel"
    ,"cooksD"
    ,"cooksDist"
    ,"OOpreds"
    ,"bb"
    ,"predLogL"
    ,"edf"
    ,"fit"
    ,"aic"
    ,"error.msg"
    ,"disty.cases"
    ,"dist1"
    ,"coef.obj"
    ,"colnammy"
    ,"offy"
    ,"mf"
    ,"type"
    ,"resids"
    ,"site.logls"
    ,"ii"
    ,"X1"
    ,"nsim"
    ,"W1"
    ,"sims"
    ,"n.sim"
    ,"pwers"
    ,"G"
    ,"inits"
    ,"outcomes"
    ,"tmp1"
    ,"tmpGrp"
    ,"tmpX"
    ,"lambda.seq"
    ,"fam"
    ,"tmp.fm"
    ,"locat.s"
    ,"my.coefs"
    ,"lastID"
    ,"tail"
    ,"df3"
    ,"tmp.fm1"
    ,"fishMod"
    ,"y"
    ,"."
    ,"MASS"
    ,"preds"
    ,"my.sd"
    ,"mult"
    ,"rcp_formula"
    ,"species_formula"
    ,"form.W"
    ,"tmp.fun"
    ,"intercepts"
    ,"form.X"
    ,"en"
    ,"root"
    ,"c1"
    ,"eta"
    ,"mu"
    ,"double.eps"
    ,"sigma1"
    ,"method"
    ,"ev"
    ,"values"
    ,"retval"
    ,"vectors"
    ,"sigsvd"
    ,"d"
    ,"v"
    ,"u"
    ,"o"
    ,"D.n"
    ,"x0"
    ,"m"
    ,"D.f0"
    ,"f"
    ,"..."
    ,"D.accur"
    ,"D.w"
    ,"D.co"
    ,"D.n.c"
    ,"eps"
    ,"macheps"
    ,"D.h"
    ,"D.deriv"
    ,"mc.fun"
    ,"D.temp.f"
    ,"jj"
    ,"D.xd"
    ,"tmp.fun.vals"
    ,"theta"
    ,"scores"
    ,"start.vals"
    ,"logl.sites"
    ,"loglOnly"
    ,"derivOnly"
    ,"RCPs"
    ,"simDat"
    ,"posts"
    ,"fm"
    ,"perms"
    ,"gtools"
    ,"classErr"
    ,"classErrRunnerUp"
    ,"postsTMP"
    ,"my.tab"
    ,"perm"
    ,"G1"
    ,"G2"
    ,"new.fm"
    ,"species"
    ,"obs.resid"
    ,"shad"
    ,"alpha.conf"
    ,"allResids"
    ,"s"
    ,"newy"
    ,"allResidsSort"
    ,"quants"
    ,"envel"
    ,"sort.resid"
    ,"empQuant"
    ,"realMeans"
    ,"realDiff"
    ,"aa"
    ,"grey"
    ,"globEnvel"
    ,"sppID"
    ,"spp.cols"
    ,"main"
    ,"fitted.scale"
    ,"loggy"
    ,"oosSizeRange"
    ,"oosDiffs"
    ,"oosWidth"
    ,"minWidth"
    ,"histy"
    ,"predlogls"
    ,"ncuts"
    ,"max.dens"
    ,"ylimmo"
    ,"breaks"
    ,"rgb"
    ,"colorRamp"
    ,"newdata"
    ,"titbit"
    ,"object2"
    ,"nboot"
    ,"my.nboot"
    ,"allCoBoot"
    ,"alphaBoot"
    ,"tauBoot"
    ,"betaBoot"
    ,"alphaIn"
    ,"tauIn"
    ,"betaIn"
    ,"gammaIn"
    ,"dispIn"
    ,"powerIn"
    ,"predCol"
    ,"ptPreds"
    ,"bootPreds"
    ,"conc"
    ,"mysd"
    ,"myContr"
    ,"boot.funny"
    ,"bootSampsToUse"
    ,"seg"
    ,"bPreds"
    ,"row.exp"
    ,"ses"
    ,"cis"
    ,"n.tot"
    ,"dat"
    ,"Call"
    ,"Distribution"
    ,"n.reorder"
    ,"all.wts"
    ,"MLstart"
    ,"my.inits"
    # ,"orderSamps"
    ,"my.fun"
    ,"dummy"
    ,"dumbOut"
    ,"capture.output"
    ,"samp.object"
    ,"flag"
    ,"tmpOldQuiet"
    ,"boot.estis"
    ,"drop.unused.levels"
    ,"stats"
    ,"optimiseDisp"
    ,"nstart"
    ,"tmpQuiet"
    ,"many.starts"
    ,"fn"
    ,"logdisp"
    ,"tmpPow"
    ,"tmpLower"
    ,"tmpUpper"
    ,"nonzero"
    ,"tmpObs"
    ,"RQR.fun"
    ,"yi"
    ,"many_yi"
    ,"F_i"
    ,"F_i_minus"
    ,"r_i"
    ,"len"
    ,"sm"
    ,"logDisps"
    ,"powers"
    ,"sppNames"
    ,"etaPi"
    ,"habis"
    ,"etaMu"
    ,"etaMu1"
    ,"etaMu2"
    ,"hh"
    ,"doPlot"
    ,"Tw.phi.func"
    ,"disp3"
    ,"spp3"
    ,"phi1"
    ,"Tw.phi.func.grad"
    ,"tmp.disp.score"
    ,"kount"
    ,"tmp.new"
    ,"tmp.old"
    ,"objective"
    ,"error"
    ,"hess"
    ,"D.accuracy"
    ,"vcov.mat"
    ,"coefMat"
    ,"summy"
    ,"emp.info"
))



