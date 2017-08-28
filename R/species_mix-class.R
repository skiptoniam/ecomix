#' @title species_mix objects
#' @rdname species_mix-classs
#' @name species_mix
#' @description Fits a finite mixture model to identify species archetype models (SAMs).
#' @details species_mix is used to fit mixtures of glms to multivariate species data. The function uses BFGS to optimise the mixture likelihood. There is the option to use EM get appropriate starting parameters. Species acts as a wrapper for fitmix.cpp that allows for easier data input. The data frames are merged into the appropriate format for the use in fitmix.cpp. Minima is found using vmmin (BFGS) and the gradients are calculated using CPPAD (auto differentiation)
#' @param formula an object of class "formula" (or an object that can be coerced to that class).
#' The response variable (left hand side of the formula) needs to be either 'presence', 'occurrence', 'abundance', 'biomass' or 'quantity' this will help specify the type of data to be modelled, if the response variable is disperate to the model distribution an error will be thrown. The dependent variables (the right hind side) of this formula specifies the dependence of the species archetype probabilities on covariates. An example formula follows something like this: cbind(spp1,spp2,spp3)~1+temperature+rainfall
#' @param model_data a List which contains named objects 'species_data': a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned. And 'covariate_data' a data frame containng the covariate data for each site. Names of columns must match that given in \code{formula}.
#' @param n_mixtures The number of mixing components (groups) to fit.
#' @param distribution The family of statistical distribution to use within the ecomix models. a  choice between "bernoulli", "poisson", "ipp" (inhomogeneous point process), "negative_binomial", "tweedie" and "gaussian" distributions are possible and applicable to specific types of data.
#' @param offset a numeric vector of length nrow(data) that is included into the model as an offset. It is included into the conditional part of the model where conditioning is performed on the SAM.
#' @param weights a numeric vector of length nrow(data) that is used as weights in the log-likelihood calculations. If NULL (default) then all weights are assumed to be identically 1. For ipp distribution - weights must be a nrow(data)*n_species matrix, which provides a species-specific background weights used to estimate the species-specific marginal likelihood.
#' @param control a list of control parameters for optimisation and calculation. See details. From \code{species_mix.control} for details on optimistaion parameters.
#' @param inits NULL a numeric vector that provides approximate starting values for species_mix coefficents. These are distribution specific, but as a minimum you'll need pis, alphas (intercepts) and betas.
#' @export
#' @examples
#' simulated_data <- simulate_species_mix_data()
#' form <- cbind(spp1,spp2,spp3) ~ 1 + x1 + x2 + x3
#' model_data <- make_mixture_data(species_data = Y, covariate_data = X)
#' fm_species_mix <- species_mix(formula,model_data=model_data,distribution='bernoulli',n_mixtures=5)

"species_mix" <- function(formula = NULL, data, n_mixtures = 3, distribution="poisson",
  offset=NULL, weights=NULL, control=species_mix.control(), standardise = FALSE){

  #the control parameters
  control <- set_control_sm(control)
  if(!control$quiet)
    message( "SAM modelling")
  call <- match.call()
  if(!is.null(formula))
    formula <- as.formula(formula)
  else{
    if(!control$quiet)
      message("There is no SAM model! Please provide a model (intercept at least) -- exitting now")
    return(NULL)
  }

  # Create model matrix
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  if(distribution=="ipp") mf$na.action <- "na.pass"
  else mf$na.action <- "na.exclude"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-1:nrow(mf)

  # get responses
  y <- model.response(mf)

  # logical matirx needed for removing NAs from response and weights.
  if(distribution=='ipp')y_is_na <- is.na(y)
  else y_is_na <- NULL
  # print(dim(y_is_na))
  # check names of reponses
  S <- check_reponse_sam(y)

  if (!S){
    if(!control$quiet)
      message("Two species have the same name -- exitting now")
    return(NULL)
  }
  if( !control$quiet)
    message( "There are ", n_mixtures, " archtypes to group the species into")

  # get model matrix
  X <- model.matrix(formula,mf)

  #get distribution
  disty.cases <- c("bernoulli","poisson","ipp","negative_binomial","tweedie","normal")
  disty <- get_distribution_sam(disty.cases, distribution)

  # get offsets and weights
  offy <- get_offset_sam(mf)
  wts <- get_weights_sam(mf,S,distribution)

  if(distribution=='ipp'){
    if(!all(colnames(y)%in%colnames(wts)))
      stop('When modelling a inhomogenous poisson point process weights colnames must match species data colnames')
    if(any(dim(y)!=dim(wts)))
      stop('Weights needs to have the same dimensions at the species data - sites x species')
    # change the response variable to be a weighted incident function 'z'.
    # y <- y/wts
  }

  s.means = NULL
  s.sds = NULL
  if (standardise == TRUE) {
    stand.X = standardise.X(X[, -1])
    X = as.matrix(cbind(1, stand.X$X))
    s.means = stand.X$dat.means
    s.sds = stand.X$dat.sds
  }

  # summarising data to console
  print_input_sam(y, X, S, formula, distribution, quiet=control$quiet)

  # used wrapper to run Piers' models.
  # fit this bad boy. bad boys, bad boys, what you gonna do when they come for you.
  if(any(distribution!=c('poisson','ipp'))){
    tmp <- fit_species_mix_wrapper(y=y, X=X, weights=wts, offset=offy, distribution=disty, G=n_mixtures, control=control, y_is_na=y_is_na, estimate_variance=control$est.var)
  } else {
    tmp <- species_mix.fit(y=y, X=X, weights=wts, offset=offy, distribution=disty, G=n_mixtures, control=control, y_is_na=y_is_na, estimate_variance=control$est.var)
    tmp$formula <- formula
    class(tmp) <- c("archetype",disty)
  }
  return(tmp)
}


#'@rdname species_mix-classs
#'@name species_mix.fit
#'@param y is a matrix genertated from \link[stats]{model.response} containing the species information. The matrix has the dimensions n_sites * n_species.
#'@param X is a design matrix of dimension n_sites * n_covariates.
#'@param G is the number of species archetypes that are being estimated.
#'@param weights is used in alternative way depending on the error distribution used. See \link[ecomix]{species_mix} for more details.
#'@param distribution the error distribution to used in species_mix estimation. Currently, 'bernoulli', 'poisson', 'ipp' (Poisson point process), 'negative_binomial' and 'tweedie' are avaliable.
#'@param offset this is a vector of site specific offsets, this might be something like area sampled at sites.
#'@param control this is a list of control parameters that alter the specifics of model fitting. See \link[ecomix]{species_mix.control} for details.
#'@param y_is_na This is a logical matrix used specifically with 'ipp' modelling - don't worry about this, it'll be worked out for you. Yay!

# test_fit_mix_ipp <- fitmix_ipp(y, X, weights, offset, G = 4, control)
# a function to fit species mix. irrespective of distribution. I should be able to wrap this around pisers' existing distributions.
"species_mix.fit" <- function(y, X, G, weights, offset, distribution, control, y_is_na=NULL){

  #use the wrapper function as a way to deal with distributions and existing fit steps.

  if(distribution == 2) tmp <- fitmix_poisson(y, X, G, weights, offset, control, y_is_na)
  if(distribution == 3) tmp <- fitmix_ipp(y, X, G, weights, offset, control, y_is_na)
  else stop('current only ipp or Poisson distribution is set up to use "species_mix.fit"')
  return(tmp)
}

#'@rdname species_mix-classs
#'@name control
#'@param quite Should any reporting be performed? Default is FALSE, for reporting.
#'@param trace int 1=model will report parameter estimates and loglikelihood at each iteration. 0=quite.
#'@param reltol function that determines the relative tolernace for model convergence. Default is quite strict.
#'@param maxit Maximum number of evaluations of the objective function allowed. Defaults to 500.
#'@param cores The number of cores to use in fitting of species mix models. These will be largely used to model the species-specific parameteres.
#'@param em.prefit Logical if TRUE the model will run a slower EM algorithim fit to find starting values.
#'@param em.steps int Default is 3, the number of EM iterations to get to starting values.
#'@param em.refit int Default is 1, number of times to refit using EM.
#'@param est.var logical if TRUE model will numerically estimate the variance covariance matrix.
#'@param residuals logical if TRUE model will estimate residuals.
#'@export

"species_mix.control" <- function(maxit = 500,
  quiet = FALSE,
  trace = 1,
  reltol = reltol_fun,
  cores = 4,
  em.prefit = TRUE,
  em.steps = 3,
  em.refit = 1,
  est.var = FALSE,
  ...){
  rval <- list(maxit = maxit, quiet = quiet, trace = trace, nreport=nreport, reltol = reltol, cores = cores,
    em.refit = em.refit, em.steps = em.steps, em.refit = em.refit, est.var = est.var)
  rval <- c(rval, list(...))
  if (is.null(rval$reltol))
    rval$reltol <- sqrt(.Machine$double.eps)
  rval
}

reltol_fun <- function(logl_n1, logl_n){
  return(abs(logl_n1 - logl_n) > (abs(logl_n1 - logl_n) / abs(logl_n)))
}

#'@rdname species_mix-classs
#'@name species_mix.predict
#'@param object is a matrix model returned from the species_mix model.
#'@param new_obs a matrix of new observations for prediction.
#'@description Predicts SAM probabilities at a series of sites. Confidence intervals can be calculated if variance-covariance matrix is estimated during species_mix model fit.
#'@examples
#'fm1 <- species_mix(form,data)
#'preds_fm1 <- predict(fm1,newdata)

"species_mix.predict" <-function (object, new_obs, ...){
  mixture.model <- object
  if (class(mixture.model)[2] == "bernoulli") {
    G <- length(mixture.model$pi)
    covar <- mixture.model$covar[-(1:(G - 1)), -(1:(G -
        1))]
    coef <- mixture.model$coef
    model.fm <- as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- model.matrix(model.fm, new_obs)
    link.fun <- make.link("logit")
    outvar <- matrix(NA, dim(X)[1], G)
    outpred <- matrix(NA, dim(X)[1], G)
    colnames(outvar) <- colnames(outpred) <- paste("G",
      1:G, sep = ".")
    for (g in 1:G) {
      lp <- as.numeric(X %*% coef[g, ])
      outpred[, g] <- link.fun$linkinv(lp)
      dhdB <- (exp(lp)/(1 + exp(lp))) * X - exp(lp)^2/((1 +
          exp(lp))^2) * X
      c2 <- covar[seq(g, dim(covar)[1], G), seq(g, dim(covar)[1],
        G)]
      for (k in 1:dim(X)[1]) {
        outvar[k, g] <- (dhdB[k, ] %*% c2) %*% (dhdB[k, ])
      }
    }
  }
  if (class(mixture.model)[2] == "negative_binomial" | class(mixture.model)[2] ==
      "tweedie") {
    G <- length(mixture.model$pi)
    covar <- mixture.model$covar[-1 * c(1:(G - 1), (dim(mixture.model$covar)[1] -
        length(mixture.model$sp.intercept) + 1):dim(mixture.model$covar)[1]),
      -1 * c(1:(G - 1), (dim(mixture.model$covar)[1] -
          length(mixture.model$sp.intercept) + 1):dim(mixture.model$covar)[1])]
    sp.int <- mixture.model$sp.intercept
    coef <- mixture.model$coef
    model.fm <- as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- cbind(model.matrix(model.fm, new_obs), 1)
    offset <- model.frame(model.fm, data = new_obs)
    offset <- model.offset(offset)
    if (is.null(offset))
      offset <- rep(0, nrow(X))
    outvar <- matrix(NA, dim(X)[1], G)
    outpred <- matrix(NA, dim(X)[1], G)
    colnames(outvar) <- colnames(outpred) <- paste("G",
      1:G, sep = ".")
    for (g in 1:G) {
      s.outvar <- matrix(NA, dim(X)[1], length(sp.int))
      s.outpred <- matrix(NA, dim(X)[1], length(sp.int))
      for (s in 1:length(sp.int)) {
        lp <- as.numeric(X %*% c(coef[g, ], sp.int[s]) +
            offset)
        s.outpred[, s] <- exp(lp)
        dhdB <- exp(lp) * X
        c2 <- covar[c(seq(g, G * (dim(X)[2] - 1), G),
          G * (dim(X)[2] - 1) + s), c(seq(g, G * (dim(X)[2] -
              1), G), G * (dim(X)[2] - 1) + s)]
        for (k in 1:dim(X)[1]) {
          s.outvar[k, s] <- (dhdB[k, ] %*% c2) %*% (dhdB[k,
            ])
        }
      }
      outpred[, g] <- apply(s.outpred * rep(mixture.model$tau[,
        g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[,
          g])
      outvar[, g] <- apply(s.outvar * rep(mixture.model$tau[,
        g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[,
          g])
    }
  }
  if (class(mixture.model)[2] == "ipp" | class(mixture.model)[2] == "poisson") {
    G <- length(mixture.model$pi)
    covar <- mixture.model$covar[-1 * c(1:(G - 1)), -1 * c(1:(G - 1))]
    sp.int <- mixture.model$sp_intercept
    coef <- mixture.model$coef
    model.fm <- as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- cbind(model.matrix(model.fm, new_obs))
    offset <- model.frame(model.fm, data = new_obs)
    offset <- model.offset(offset)
    if (is.null(offset))
      offset <- rep(0, nrow(X))
    outvar <- matrix(NA, dim(X)[1], G)
    outpred <- matrix(NA, dim(X)[1], G)
    colnames(outvar) <- colnames(outpred) <- paste("G", 1:G, sep = ".")
    for (g in 1:G) {
      s.outvar <- matrix(NA, dim(X)[1], length(sp.int))
      s.outpred <- matrix(NA, dim(X)[1], length(sp.int))
      for (s in 1:length(sp.int)) {
        lp <- as.numeric(X %*% c(sp.int[s],coef[g, ]) + offset)
        s.outpred[, s] <- exp(lp)
        dhdB <- exp(lp) * X
        ## pull out the sp intercept and mixing component covariates.
        c2 <- covar[c(G * (dim(X)[2] - 1) + s,seq(g, G * (dim(X)[2] - 1), G)), c((G * (dim(X)[2] - 1) + s),seq(g, G * (dim(X)[2] - 1), G))]
        for (k in 1:dim(X)[1]) {
          s.outvar[k, s] <- (dhdB[k, ] %*% c2) %*% (dhdB[k, ])
        }
      }
      outpred[, g] <- apply(s.outpred * rep(mixture.model$tau[, g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[, g])
      outvar[, g] <- apply(s.outvar * rep(mixture.model$tau[, g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[, g])
    }
  }
  list(fit = outpred, se.fit = sqrt(outvar))
}

"print_input_sam" <- function(y, X, S, formula, distribution, quiet=FALSE){
  if( quiet)
    return( NULL)
  n.tot <- nrow(y)
  if(distribution=='ipp'){
    n_pres <- sum(unlist(y)==1,na.rm=TRUE)
    n_bkgrd <- sum(unlist(y[,1])==0,na.rm=TRUE)
    message("There are ", n_pres, " presence observations for ", S," species")
    message("There are ", n_bk_pts, " background (integration) points for each of the ", S," species")
  } else {
    message("There are ", nrow(X), " site observations for ", S," species")
  }
  formula[[2]] <- NULL
  message("The model for the SAM is ", Reduce( "paste", deparse( formula)))
  message("You are implementing a ", distribution, " SAM.")
}

"get_distribution_sam" <- function( disty.cases, dist1) {
  error.msg <- paste( c( "Distribution not implemented. Options are: ", disty.cases, "-- Exitting Now"), collapse=" ")
  disty <- switch( dist1, "bernoulli" = 1,"poisson" = 2,"ipp"=3,"negative_binomial" = 4,"tweedie" = 5,"gaussian" = 6,{stop( error.msg)} )
  return( disty)
}

"get_X_sam" <- function(form.SAM, mf.X){
  form.X <- form.SAM
  form.X[[2]] <- NULL
  form.X <- as.formula(form.X)
  X <- model.matrix(form.X, mf.X)
  return( X)
}

"check_reponse_sam" <-function(outs) {
  nam <- colnames( outs)
  if( length( nam) == length( unique( nam)))
    return( length( nam))
  else
    return( FALSE)
}
