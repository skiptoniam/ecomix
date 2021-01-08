##### Main species mix functions to export #####
#' @title This is how you fit a species archetype model (SAM) in ecomix.
#' @rdname species_mix
#' @name species_mix
#' @description Fits a mixture-of-regressions to identify species archetype
#' models (SAMs).
#' @details species_mix is used to fit mixtures of glms to multivariate
#' species data. The function uses BFGS to optimize the mixture likelihood.
#' There is the option to use ECM algorithm to get appropriate starting
#' parameters. `species_mix` acts as a wrapper for species_mix.fit
#' that allows for easier data input. The data frames are merged into
#' the appropriate format for the use in species_mix.fit.
#' Minima is found using vmmin (BFGS). Currently 'bernoulli', 'binomial',
#' 'poisson', 'ippm' (inhomogeneous Poisson point process), 'negative.binomial',
#' 'tweedie' and 'gaussian' distributions can be fitted using the species_mix
#' function.
#' @param archetype_formula an object of class "formula" (or an object that can
#' be coerced to that class). The response variable (left hand side of the
#' formula) needs to be either 'occurrence', 'abundance',
#' 'biomass' or 'quantity' data. The type of response data will help specify
#' the type of error distribution to be used. The dependent variables
#' (the right hind side) of this formula specifies the dependence of the
#' species archetype probabilities on covariates. For all model the basic
#' formula structure follows something like this:
#' cbind(spp1,spp2,spp3)~1+temperature+rainfall
#' @param species_formula an object of class "formula" (or an object that can be
#' coerced to that class). The right hand side of this formula specifies the
#' dependence of the species"'" data on covariates (typically different
#' covariates to \code{archetype_formula} to avoid confusing confounding).
#' Current the formula is set at ~ 1 by default for species-specific intercepts
#' for the archetype models. If you include a species specific formula which has
#' more than an intercept you will be fitting a partial species archetype model
#' which has species specific covariates and archetype specific covariates.
#' @param all_formula an object of class "formula", which is meant to represent
#'  a constant single set of covariates across all species and groups, typically
#'  you might use this an alternative to an offset, where there might be some
#'  bias in the data which is relatively constant and might arise as an
#'  artefact of how the data was collected.
#' @param data a matrix or data.frame which contains the 'species_data'
#' matrix, a const and the covariates in the structure of spp1, spp2, spp3,
#' const, temperature, rainfall. dims of matrix should be
#' nsites*(nspecies+const+covariates).
#' @param nArchetypes The number of archetypes (mixing components/groups) to estimate from the data.
#' @param family The family of statistical family to use within
#' the ecomix models. a  choice between "bernoulli", "poisson", "ippm",
#' "negative.binomial" and "gaussian" families are possible and applicable
#' to specific types of data.
#' @param offset a numeric vector of length nrow(data) (n sites) that is
#' included into the model as an offset. It is included into the conditional
#' part of the model where conditioning is performed on the SAM.
#' @param weights a numeric vector of length ncol(Y) (n species) that is used
#' as weights in the log-likelihood calculations. If NULL (default) then all
#' weights are assumed to be identically 1. Because we are estimating the
#' log-likelihood over species (rather than sites), the weights should be a
#' vector n species long.
#' @param bb_weights a numeric vector of n species long. This is used for
#' undertaking a Bayesian Bootstrap. See 'vcov.species_mix' for more details.
#' @param size The size of the sample for a binomial model (defaults to 1).
#' @param power The power parameter for a Tweedie model. Default is 1.6, and this is assigned to all species
#' @param control a list of control parameters for optimisation and calculation.
#' See details.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficients. These are family specific, but at a
#' minimum you will need pis (additive_logitic transformed), alpha
#' (intercepts) and beta (mixing coefs).
#' @param standardise Boolean. If TRUE, standardise the covariate data.
#' @param titbits either a boolean or a vector of characters. If TRUE
#' (default for species_mix(qv)), then some objects used in the estimation of
#' the model"'"s parameters are returned in a list entitled "titbits" in the
#' model object. Some functions, for example plot.species_mix(qv) and
#' predict.species_mix(qv), will require some or all of these pieces of
#' information. If titbits=FALSE (default for species_mix.multifit(qv)),
#' then an empty list is returned. If a character vector, then just those
#' objects are returned. Possible values are:"Y" for the outcome matrix, "X"
#' for the model matrix for the SAM model, "offset" for the offset in the model,
#' "site_spp_weights" for the model weights, "archetype_formula" for the formula
#' for the SAMs, "species_formula" for the formula for the species-specific
#' model, "control" for the control arguments used in model fitting, "family" for
#' the conditional distribution of the species data. Care needs to be taken when
#' using titbits=TRUE in species_mix.multifit(qv) calls as titbits is created
#' for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then
#' setting titbits=TRUE may give users problems with memory.
#' @importFrom graphics abline hist legend lines matplot par plot points polygon
#'  rect
#' @importFrom stats as.formula binomial cooks.distance cov cutree dbinom family dnbinom dnorm dpois make.link coef glm.fit fitted gaussian glm hclust lm logLik model.matrix model.frame model.offset model.response model.weights pbinom pnbinom pnorm poisson ppois predict qnorm qqnorm quantile rbinom residuals rgamma rnbinom rnorm rpois runif sd uniroot update update.formula nlminb optimise
#' @export
#' @examples
#' \donttest{
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
#'  }

"species_mix" <- function(archetype_formula = NULL,
                          species_formula = stats::as.formula(~1),
                          all_formula = NULL,
                          data, nArchetypes = 3,
                          family="bernoulli", offset=NULL,
                          weights=NULL, bb_weights=NULL, size = NULL, power=1.6,
                          control=list(), inits=NULL, standardise = FALSE, titbits = TRUE){

  data <- as.data.frame(data)
  control <- set_control_sam(control)
  if(!control$quiet)
    message( "SAM modelling")
  call <- match.call()
  if(!is.null(archetype_formula)){
    archetype_formula <- stats::as.formula(archetype_formula)
  } else{
    if(!control$quiet)
      message("There is no SAM model!\n
              Please provide an archetype_formula -- exitting now")
    return(NULL)
  }

  which_mix <- check_species_formula(species_formula)
  if(!is.null(species_formula))
    species_formula <- stats::as.formula(species_formula)

    mf <- match.call(expand.dots = FALSE)
  if(family=="ippm"){
    m <- match(c("data","offset"), names(mf), 0L)
  } else {
    m <- match(c("data","offset","weights"), names(mf), 0L)
  }

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  if(family=="ippm"){
    mf$na.action <- "na.pass"
  } else {
    mf$na.action <- "na.exclude"
  }

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))
  dat <- clean_data_sam(mf, archetype_formula, species_formula, all_formula, family)

  # get responses
  y <- stats::model.response(dat$mf.X)

  # logical matirx needed for removing NAs from response and weights.
  y_is_na <- is.na(y)

  # check names of reponses
  S <- check_reponse_sam(y)

  if (!S){
    if(!control$quiet)
      message("Two species have the same name -- exitting now")
    return(NULL)
  }
  if( !control$quiet)
    message( "There are ", nArchetypes, " archetypes to group the species into")

  # get archetype model matrix
  X <- get_X_sam(archetype_formula = archetype_formula, mf.X = dat$mf.X)

  # what is the W matrix (species covariates)
  W <- get_W_sam(species_formula = species_formula, mf.W = dat$mf.W)

  # what is the U matrix (species covariates)
  U <- get_U_sam(all_formula = all_formula, mf.U = dat$mf.U)

  # standarise data if requested.
  x.means <- NULL
  x.sds <- NULL
  w.means <- NULL
  w.sds <- NULL
  if (standardise == TRUE) {
    stand.X <- standardise.X(X)
    X <- as.matrix(stand.X$X)
    x.means <- stand.X$dat.means
    x.sds <- stand.X$dat.sds
    if(ncol(W)>1){
      stand.W <- standardise.W(W[,-1,drop=FALSE])
      W <- as.matrix(cbind(1,stand.W$W))
      w.means <- stand.W$dat.means
      w.sds <- stand.W$dat.sds
    }
    if(!is.null(U)){
      stand.U <- standardise.U(U)
      U <- as.matrix(stand.U$U)
      u.means <- stand.U$dat.means
      u.sds <- stand.U$dat.sds
    }
  }

  #get family
  disty_cases <- c("bernoulli","poisson","ippm",
                   "negative.binomial","tweedie",
                   "gaussian","binomial")

  disty <- get_family_sam(disty_cases, family)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names,
                                               family)
  spp_weights <- check_spp_weights(bb_weights,S)

  size <- check_size_binomial(size,nrow(dat$mf.X))

  ## check powers
  powers <- get_power_sam(disty,power,S)

  if(family=='ippm'){
    if(!all(colnames(y)==colnames(site_spp_weights))){
      stop(message('When modelling a inhomogeneous poisson point process model,
               \n species data colnames must match weights colnames.\n\n
               Species data colnames from "data" are:\n',colnames(y),'.\n\n
               While the colnames of the weights are:\n',
               colnames(site_spp_weights),'\n'))
    }
    if(any(dim(y)!=dim(site_spp_weights))){
      stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
    }
  }

  # summarising data to console
  if(!control$quiet) print_input_sam(y, X, W, U, S, archetype_formula,
                                     species_formula, all_formula,
                                     family,
                                     quiet=control$quiet)

  # fit species mix.
  G <- nArchetypes
  tmp <- species_mix.fit(y=y, X=X, W=W, U=U, G=G, S=S,
                         spp_weights=spp_weights,
                         site_spp_weights=site_spp_weights,
                         offset=offset, disty=disty, y_is_na=y_is_na,
                         size=size, powers=powers, control=control, inits=inits)

  tmp$family <- disty_cases[disty]

  tmp$S <- S;
  tmp$G <- G;
  tmp$npx <- ncol(X);
  tmp$npw <- ifelse(ncol(W)>1,ncol(W),0);
  tmp$npu <- ifelse(!is.null(U),ncol(U),0);

  if(nArchetypes==1){
    tmp$pis <- tmp$pis
  }else{
    tmp$pis <- additive_logistic(tmp$eta)
  }

  #get logls from parameters
  #calc posterior porbs and pis.
  if(nArchetypes>1){
    fits <- tmp$coefs
    logls_mus <- get_logls_sam(y = y, X = X, W = W, U = U, G = G,
                               S = S, spp_weights = spp_weights,
                               site_spp_weights = site_spp_weights,
                               offset = offset, y_is_na = y_is_na,
                               disty = disty, size = size,
                               powers = powers, control=control,
                               fits = fits, get_fitted = FALSE)
    tmp$taus <- get_taus(tmp$pis,logls_mus$logl_sp,G,S)
    tmp$pis <- colSums(tmp$taus)/S
  }

  #Information criteria
  tmp <- calc_info_crit_sam(tmp)

  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_sam(titbits, y, X, W, U, spp_weights,
                                 site_spp_weights, offset, y_is_na, size, powers,
                                 archetype_formula, species_formula, all_formula,
                                 control, disty_cases[disty])

  # remove large annoying object if titbits == FALSE
  if(!titbits)
    tmp <- titbit_cleanup(tmp)

  class(tmp) <- c("species_mix")
  return(tmp)
}

#'@title species_mix.fit
#'@rdname species_mix.fit
#'@name species_mix.fit
#' @description species_mix.fit is similar to glm.fit and does all the heavy lifting when it
#' comes to estimating species mix models. If you are unfamilar with how to use glm.fit it is
#' recommended that you use species_mix which is the user friendly wrapper around this
#' function.
#'@param y is a matrix generated from model.response containing the species information. The matrix has the dimensions n_sites * n_species.
#'@param X is a design matrix for the archetype_formula dimension n_sites * n_covariates.
#'@param W is a design matrix for species_formula and will be implemented if species_formula has covariates.
#'@param U is a design matrix for all_formula and will be implemented if not NULL.
#'@param G is the number of species archetypes that are being estimated.
#'@param S is the number of species to be modelled (this will be calculated internally in species_mix())
#'@param spp_weights These are weights on the species logls and are specifically used in the Bayesian Bootstrap.
#'@param site_spp_weights These are site and species specific weights. For most distributions these will be the same across all species. But this form is required to correctly estiamte the IPPMs. See \link[ecomix]{species_mix} for more details.
#'@param offset this is a vector of site specific offsets, this might be something like area sampled at sites.
#'@param y_is_na This is a logical matrix used specifically with 'ippm' modelling - don't worry about this, it'll be worked out for you. Yay!
#'@param disty the error distribution to used in species_mix estimation. Currently, 'bernoulli', 'poisson', 'ippm' (Poisson point process), 'negative.binomial' and 'guassian' are available - internal conversion of distribution to a integer.
#'@param size The size of each of binomial sample at each site. Length should be the number of sites.
#'@param powers The power parameters for the Tweedie distribution.
#'@param control this is a list of control parameters that alter the specifics of model fitting.
#'@param inits This will be a vector of starting values for species_mix (i.e you've fitted a model and want to refit it).
#'@export

"species_mix.fit" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
                              offset, y_is_na, disty, size, powers, control, inits=NULL){

  if(G==1){
    tmp <- fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights,
                       offset, y_is_na, G, S, disty, size, powers,
                       control=control)
    return(tmp)
  }

  if(is.null(inits)){
    starting_values  <- starting_values_wrapper(y = y, X = X, W = W, U = U,
                                                site_spp_weights = site_spp_weights,
                                                offset = offset,
                                                y_is_na = y_is_na,
                                                G = G, S = S,
                                                disty = disty,
                                                size = size,
                                                powers = powers,
                                                control = control)

  } else {
    if(!control$quiet)message('Be careful! You are using your own initial starting values to optimise the species_mix model.')
    inits <- setup_inits_sam(inits, S=S, G=G, X=X, W=W, U=U, disty, return_list = TRUE)
    print(inits)
    starting_values <- inits
  }

  tmp <- sam_optimise(y = y, X = X, W = W, U = U, offset = offset,
                      spp_weights = spp_weights,
                      site_spp_weights = site_spp_weights,
                      y_is_na = y_is_na, S = S, G = G, disty = disty,
                      size = size, powers = powers, start_vals = starting_values,
                      control = control)

  return(tmp)
}

#'@rdname species_mix.multifit
#'@name species_mix.multifit
#'@title species_mix.multifit
#'@description This version of species mix is useful for fitting models which have complex likelihoods.
#'The multiple starts will enable optimisation of the loglikelihood using multiple starts.
#' @param archetype_formula an object of class "formula" (or an object that can
#' be coerced to that class). The response variable (left hand side of the
#' formula) needs to be either 'occurrence', 'abundance',
#' 'biomass' or 'quantity' data. The type of reponse data will help specify
#' the type of error distribution to be used. The dependent variables
#' (the right hind side) of this formula specifies the dependence of the
#' species archetype probabilities on covariates. For all model the basic
#' formula structure follows something like this:
#' cbind(spp1,spp2,spp3)~1+temperature+rainfall
#' @param species_formula an object of class "formula" (or an object that can be
#' coerced to that class). The right hand side of this formula specifies the
#' dependence of the species"'" data on covariates (typically different
#' covariates to \code{archetype_formula} to avoid confusing confounding).
#' Current the formula is set at ~ 1 by default for species-specific intercepts
#' for the archetype models. If you include a species specific formula which has
#' more than an intercept you will be fitting a partial species archetype model
#' which has species specific covariates and archetype specific covariates.
#' @param all_formula an object of class "formula", which is meant to represent
#'  a constant single set of covariates across all species and groups, typically you might
#'  use this an alternative to an offset, where there might be some bias in the
#'  data which is relatively constant and might arise as an artefact of how the data was collected.
#' @param data a matrix of dataframe which contains the 'species_data'
#' matrix, a const and the covariates in the strucute of spp1, spp2, spp3,
#' const, temperature, rainfall. dims of matirx should be
#' nsites*(nspecies+const+covariates).
#' @param nArchetypes The number of mixing components (groups) to fit.
#' @param family The family of statistical distribution to use within
#' the ecomix models. a  choice between "bernoulli", "poisson", "ippm",
#' "negative.binomial" and "gaussian" distributions are possible and applicable
#' to specific types of data.
#' @param offset a numeric vector of length nrow(data) (n sites) that is
#' included into the model as an offset. It is included into the conditional
#' part of the model where conditioning is performed on the SAM.
#' @param weights a numeric vector of length ncol(Y) (n species) that is used
#' as weights in the log-likelihood calculations. If NULL (default) then all
#' weights are assumed to be identically 1. Because we are estimating the
#' log-likelihood over species (rather than sites), the weights should be a
#' vector n species long.
#' @param bb_weights a numeric vector of n species long. This is used for
#' undertaking a Bayesian Bootstrap. See 'vcov.species_mix' for more details.
#' @param size The size of the sample for a binomial model (defaults to 1).
#' @param power The power parameter for a Tweedie model. Default is 1.6, and
#' this is assigned to all species
#' @param control a list of control parameters for optimisation and calculation.
#' See details.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficients. These are distribution specific, but at a
#' minimum you will need pis (additive_logitic transformed), alpha
#' (intercepts) and beta (mixing coefs).
#' @param standardise Boolean. If TRUE, standardise the covariate data.
#' @param titbits either a boolean or a vector of characters. If TRUE
#' (default for species_mix(qv)), then some objects used in the estimation of
#' the model"'"s parameters are returned in a list entitled "titbits" in the
#' model object. Some functions, for example plot.species_mix(qv) and
#' predict.species_mix(qv), will require some or all of these pieces of
#' information. If titbits=FALSE (default for species_mix.multifit(qv)),
#' then an empty list is returned. If a character vector, then just those
#' objects are returned. Possible values are:"Y" for the outcome matrix, "X"
#' for the model matrix for the SAM model, "offset" for the offset in the model,
#' "site_spp_weights" for the model weights, "archetype_formula" for the formula
#' for the SAMs, "species_formula" for the formula for the species-specific
#' model, "control" for the control arguments used in model fitting, "family" for
#' the conditional distribution of the species data. Care needs to be taken when
#' using titbits=TRUE in species_mix.multifit(qv) calls as titbits is created
#' for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then
#' setting titbits=TRUE may give users problems with memory.
#'@param nstart for species_mix.multifit only. The number of random starts to
#'perform for re-fitting. Default is 10, which will need increasing for serious use.
#'@param mc.cores for species_mix.multifit only. The number of cores to spread
#' the re-fitting over.
#'@export
#'@examples
#' \donttest{
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
#' fmods <- species_mix.multifit(archetype_formula = sam_form,
#' species_formula = sp_form, data=simulated_data, family = 'bernoulli',
#'  nstart = 10, nArchetypes=3)
#' }
"species_mix.multifit" <- function(archetype_formula = NULL,
                                   species_formula = stats::as.formula(~1),
                                   all_formula = NULL,
                                   data, nArchetypes = 3,
                                   family="bernoulli", offset=NULL,
                                   weights=NULL, bb_weights=NULL, size = NULL,
                                   power=1.6,
                                   control=list(), inits=NULL,
                                   standardise = FALSE, titbits = TRUE,
                                   nstart = 10,
                                   mc.cores = 1){

  data <- as.data.frame(data)
  control <- set_control_sam(control)
  if(!control$quiet)
    message( "SAM modelling")
  call <- match.call()
  if(!is.null(archetype_formula)){
    archetype_formula <- stats::as.formula(archetype_formula)
  } else{
    if(!control$quiet)
      message("There is no SAM model!\n
              Please provide an archetype_formula -- exitting now")
    return(NULL)
  }

  which_mix <- check_species_formula(species_formula)
  if(!is.null(species_formula))
    species_formula <- stats::as.formula(species_formula)

  mf <- match.call(expand.dots = FALSE)
  if(family=="ippm"){
    m <- match(c("data","offset"), names(mf), 0L)
  } else {
    m <- match(c("data","offset","weights"), names(mf), 0L)
  }

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  if(family=="ippm"){
    mf$na.action <- "na.pass"
  } else {
    mf$na.action <- "na.exclude"
  }

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))
  dat <- clean_data_sam(mf, archetype_formula, species_formula, all_formula, family)

  # get responses
  y <- stats::model.response(dat$mf.X)

  # logical matirx needed for removing NAs from response and weights.
  y_is_na <- is.na(y)

  # check names of reponses
  S <- check_reponse_sam(y)

  if (!S){
    if(!control$quiet)
      message("Two species have the same name -- exitting now")
    return(NULL)
  }
  if( !control$quiet)
    message( "There are ", nArchetypes, " archetypes to group the species into")

  # get archetype model matrix
  X <- get_X_sam(archetype_formula = archetype_formula, mf.X = dat$mf.X)

  # what is the W matrix (species covariates)
  W <- get_W_sam(species_formula = species_formula, mf.W = dat$mf.W)

  # what is the U matrix (species covariates)
  U <- get_U_sam(all_formula = all_formula, mf.U = dat$mf.U)

  # standarise data if requested.
  x.means <- NULL
  x.sds <- NULL
  w.means <- NULL
  w.sds <- NULL
  if (standardise == TRUE) {
    stand.X <- standardise.X(X)
    X <- as.matrix(stand.X$X)
    x.means <- stand.X$dat.means
    x.sds <- stand.X$dat.sds
    if(ncol(W)>1){
      stand.W <- standardise.W(W[,-1,drop=FALSE])
      W <- as.matrix(cbind(1,stand.W$W))
      w.means <- stand.W$dat.means
      w.sds <- stand.W$dat.sds
    }
    if(!is.null(U)){
      stand.W <- standardise.U(U)
      U <- as.matrix(stand.U$U)
      u.means <- stand.U$dat.means
      u.sds <- stand.U$dat.sds
    }
  }

  #get family
  disty_cases <- c("bernoulli","poisson","ippm",
                   "negative.binomial","tweedie",
                   "gaussian","binomial")

  disty <- get_family_sam(disty_cases, family)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names,
                                               family)
  spp_weights <- check_spp_weights(bb_weights,S)

  size <- check_size_binomial(size,nrow(dat$mf.X))

  ## check powers
  powers <- get_power_sam(disty,power,S)

  if(family=='ippm'){
    if(!all(colnames(y)==colnames(site_spp_weights))){
      stop(message('When modelling a inhomogeneous poisson point process model,
               \n species data colnames must match weights colnames.\n\n
               Species data colnames from "data" are:\n',colnames(y),'.\n\n
               While the colnames of the weights are:\n',
                   colnames(site_spp_weights),'\n'))
    }
    if(any(dim(y)!=dim(site_spp_weights))){
      stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
    }
  }

  # summarising data to console
  if(!control$quiet) print_input_sam(y, X, W, U, S, archetype_formula,
                                     species_formula, all_formula,
                                     family,
                                     quiet=control$quiet)

  tmp_fun <- function(x){
    if( !control$quiet & nstart>1)
      setTxtProgressBar(pb, x)
        tmpQuiet <- control$quiet
        control$quiet <- TRUE
      # fit species mix.
        G <- nArchetypes
        tmp <- species_mix.fit(y=y, X=X, W=W, U=U, G=G, S=S,
                               spp_weights=spp_weights,
                               site_spp_weights=site_spp_weights,
                               offset=offset, disty=disty, y_is_na=y_is_na,
                               size=size, powers=powers, control=control, inits=inits)

        tmp$family <- disty_cases[disty]

        if(nArchetypes==1){
          tmp$pis <- tmp$pis
        }else{
          tmp$pis <- additive_logistic(tmp$eta)
        }

        # get logls from parameters

        #calc posterior porbs and pis.
        if(nArchetypes>1){
          fits <- tmp$coefs
          logls_mus <- get_logls_sam(y = y, X = X, W = W, U = U, G = G,
                                     S = S, spp_weights = spp_weights,
                                     site_spp_weights = site_spp_weights,
                                     offset = offset, y_is_na = y_is_na,
                                     disty = disty, size = size,
                                     powers = powers, control=control,
                                     fits = fits, get_fitted = FALSE)
          tmp$taus <- get_taus(tmp$pis,logls_mus$logl_sp,G,S)
          tmp$pis <- colSums(tmp$taus)/S
        }

        #Information criteria
        tmp <- calc_info_crit_sam(tmp)

        #titbits object, if wanted/needed.
        tmp$titbits <- get_titbits_sam(titbits, y, X, W, U, spp_weights,
                                       site_spp_weights, offset, y_is_na, size, powers,
                                       archetype_formula, species_formula, all_formula,
                                       control, disty_cases[disty])

        # remove large annoying object if titbits == FALSE
        if(!titbits)
          tmp <- titbit_cleanup(tmp)

        class(tmp) <- c("species_mix")
      return( tmp)
   }

  #    require( parallel)
  if( !control$quiet & nstart>1)
    pb <- txtProgressBar(min = 1, max = nstart, style = 3, char = "c[_] ")

   #Fit the model many times
   many_starts <- plapply(seq_len(nstart), tmp_fun,
                         .parallel = mc.cores,
                          .verbose = !control$quiet)

   class(many_starts) <- c("species_mix.multifit")

   return(many_starts)
}

#' @rdname species_mix.bootstrap
#' @name species_mix.bootstrap
#' @title species_mix.bootstrap
#' @param object A species_mix model object
#' @param nboot The number of boostraps to run.
#' @param type the type of bootstrap to use, options are  "SimpleBoot" which is  a parameteric bootstrap, or "BayesBoot"
#' @param mc.cores The number of core to fit. The default is 1.
#' @param quiet If TRUE, do not print progress of bootstrap.
#' @importFrom stats vcov
#' @export

"species_mix.bootstrap" <-function (object, nboot=1000, type="BayesBoot",
                                    mc.cores=1, quiet=FALSE){

  if(object$titbits$family=='ippm'){
    message('Using parameteric bootstrap to family of parameters for IPPM model\n')
    if(is.null(object$vcov)){
      message('Variance-Covariance matrix is missing, will attempt to calculate now, this will be slow...')
      object$vcov <- vcov(object)
    }

    ## this will do a parameteric bootstrap and return a matrix boots*params
    boot.estis <- mvtnorm::rmvnorm(nboot,unlist(object$coefs),object$vcov)
    tmpOldQuiet <- object$titbits$control$quiet

  } else {

  if (nboot < 1)
    stop( "No Bootstrap samples requested.  Please set nboot to something > 1.")
  if( ! type %in% c("BayesBoot","SimpleBoot"))
    stop( "Unknown Bootstrap type, choices are BayesBoot and SimpleBoot.")
  n.reorder <- 0
  object$titbits$control$optimise <- TRUE #just in case it was turned off
  if(object$titbits$family=='ippm')
    stop('IPPM vcov matrix needs to estimated using FiniteDifference method.\n')

  if( !quiet){
    chars <- c("><(('> ","_@_'' ")
    pb <- txtProgressBar(min = 1, max = nboot, style = 3, char = chars[sample(length(chars),1)])
  }

  if( type == "SimpleBoot"){
    all.wts <- matrix( sample( 1:object$S, nboot*object$S, replace=TRUE), nrow=nboot, ncol=object$S)
    tmp <- apply( all.wts, 1, table)
    all.wts <- matrix( 0, nrow=nboot, ncol=object$S)
    for( ii in seq_along( tmp))
      all.wts[ii, as.numeric( names( tmp[[ii]]))] <- tmp[[ii]]
  }
  if( type == "BayesBoot")
    all.wts <- object$S * gtools::rdirichlet( nboot, rep( 1, object$S))
  my.inits <- object$coef

  tmpOldQuiet <- object$titbits$control$quiet
  object$titbits$control$quiet <- TRUE

  my.fun <- function(dummy){
    if( !quiet) setTxtProgressBar(pb, dummy)
    disty_cases <- c("bernoulli","binomial","poisson", "ippm", "negative.binomial", "tweedie", "gaussian")
    disty <- get_family_sam(disty_cases, object$family)
    dumbOut <- capture.output(
      samp.object <- species_mix.fit(y=object$titbits$Y,
                                     X=object$titbits$X,
                                     W=object$titbits$W,
                                     U=object$titbits$U,
                                     offset = object$titbits$offset,
                                     spp_weights = all.wts[dummy,,drop=TRUE],
                                     site_spp_weights = object$titbits$site_spp_weights,
                                     G = object$G,
                                     S = object$S,
                                     y_is_na = object$titbits$y_is_na,
                                     disty = disty,
                                     size = object$titbits$size,
                                     powers = object$titbits$powers,
                                     control = object$titbits$control,
                                     inits = my.inits))
    return( unlist( samp.object$coef))
  }

  flag <- TRUE
  if( Sys.info()['sysname'] == "Windows" | mc.cores==1){
    boot.estis <- matrix(NA, nrow = nboot, ncol = length(unlist(object$coef)))
    for (ii in 1:nboot) {
      if( !quiet)
        setTxtProgressBar(pb, ii)
      boot.estis[ii, ] <- my.fun( ii)
    }
    flag <- FALSE
  }

  if( flag){
    tmp <- plapply(seq_len(nboot), my.fun, .parallel = mc.cores)
    boot.estis <- do.call( "rbind", tmp)
  }
  }
  object$titbits$control$quiet <- tmpOldQuiet
  class( boot.estis) <- "species_mix.bootstrap"
  return( boot.estis)
}

#' @title Simulate species mix data for model fitting.
#' @rdname species_mix.simulate
#' @name species_mix.simulate
#' @param archetype_formula formula to simulate species_mix data, needs to have
#' the format: cbind(spp1,spp2,spp3,...,sppN)~1 + x1 + x2
#' @param species_formula formula to simulate species_mix species-specific
#' responses, e.g: ~1
#' @param all_formula formula to simulate biases in the data
#' @param data a matrix of variables to simulate data from.
#' @param nArchetypes number of groups to simulate.
#' @param alpha coefficients for each species archetype. vector S long.
#' @param beta coefficients for each species archetype. Matrix of G x number of
#'  parameters. Each row is a different species archetype.
#' @param gamma coefficients for each species archetype. Matrix of S x number of
#'  parameters. Each row is a different species archetype.
#' @param delta coefficients for all_formula, these should describe overall
#'  biases in the dataset.
#' @param logTheta coefficients for the dispersion variables for negative.binomial
#' and gaussian distributions - should be number of species long and on the
#' natural log scale.
#' @param size Is for the binomial model and this represents the number of binomial trials per site, can be fixed or vary.
#' @param powers Is the power parameter for Tweedie distribution.
#' @param family Which statistical distribution to simulate data for.
#'  'bernoulli','binomial', 'gaussian', 'ippm', 'negative.binomial' and 'poisson'.
#' @param offset used to offset sampling effort for abundance data (log link function).
#' @export
#' @examples
#' \donttest{
#' archetype_formula <- stats::as.formula(paste0('cbind(',paste(paste0('spp',
#' 1:20),collapse = ','),")~1+x1+x2"))
#' species_formula <- stats::as.formula(~1)
#' beta <- matrix(c(-3.6,0.5,
#'                  -0.9,1.0,
#'                   0.9,-2.9,
#'                   2.2,5.4),
#'                 4,2,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),
#'                   x1=stats::runif(100,0,2.5),
#'                   x2=stats::rnorm(100,0,2.5))
#' simulated_data <- species_mix.simulate(archetype_formula,species_formula,
#'                                        data=dat, nArchetypes = 4, beta=beta,
#'                                        family="bernoulli")
#' }

## need to update this to take the new formula framework and simulate ippm data.
"species_mix.simulate" <-  function(archetype_formula,
                                    species_formula,
                                    all_formula=NULL,
                                    data,
                                    offset = NULL,
                                    nArchetypes = 3,
                                    alpha=NULL,
                                    beta=NULL,
                                    gamma=NULL,
                                    delta=NULL,
                                    logTheta=NULL,
                                    powers=NULL,
                                    size=NULL,
                                    family = "bernoulli"){

  #update the formula to old format.
  if(!is.null(archetype_formula))
    archetype_formula <- stats::as.formula(archetype_formula)
  if(!is.null(species_formula))
    species_formula <- stats::as.formula(species_formula)

  sam_org <- archetype_formula
  spp_org <- species_formula

  if(!is.null(all_formula)){
    all_formula <- stats::as.formula(all_formula)
    all_org <- all_formula

    #drop intercept from all
    all_formula <- update(all_formula,~.-1)
    U <- stats::model.matrix(all_formula, data)
    npu <- ncol(U)
  } else {
    all_org <- NULL
    npu <- 0
    U <- NULL
  }


  # how many species to simulate???
  S <- length(archetype_formula[[2]])-1
  G <- nArchetypes

  archetype_formula <- update(archetype_formula,y~.-1)
  archetype_formula[[2]] <- NULL

  X <- stats::model.matrix(archetype_formula, data)
  W <- stats::model.matrix(species_formula, data)
  if(is.null(offset)) offset <- rep(0,nrow(X))

  n <- nrow(X)
  npx <- ncol(X)
  npw <- ncol(W) - 1

  size  <- check_size_binomial(size, n)

  if (is.null(alpha)) {
    message("Random alpha from normal (-1,0.5) distribution")
    alpha <- rnorm(S,-1,0.5)
  }
  if (is.null(beta) | length(beta) != (npx) * G) {
    message("Random values for beta")
    beta <- rnorm(npx*(G))
  }
  beta <- matrix(as.numeric(beta), nrow = G)
  if(( is.null(gamma) | length( gamma) != S * npw)){
    if( npw != 0){
      message("Random values for gamma")
      gamma <- rnorm( npw*S)
      gamma <- matrix( as.numeric(gamma), nrow=S, ncol=npw)
    } else {
      gamma <- NULL
    }
  } else {
    gamma <- matrix( as.numeric( gamma), nrow=S)
  }

  if(!npu==0){
    if(is.null(delta) | length(delta) != (npu)){
      message("Random values for delta")
      delta <- rnorm(npu)
    }
  }

  if( family == "negative.binomial" & (is.null(logTheta) | length(logTheta) != S)){

    message( "Random values for overdispersions")
    logTheta <- log( 1 + rgamma( n=S, shape=1, scale=0.75))
  }
  if (family == "tweedie" & (is.null(logTheta) | length(logTheta) != S)) {
    message("Random values for species' dispersion parameters")
    logTheta <- log(1 + rgamma(n = S, shape = 1, scale = 0.75))
  }
  if (family == "tweedie" & (is.null(powers) | length(powers) != S)) {
    message("Power parameter assigned to 1.6 for each species")
    powers <- rep(1.6, S)
  }
  if( family=="gaussian"){
    if((is.null( logTheta) | length( logTheta) != S)){
    message( "Random values for species' variance parameters")
    logTheta <- log( 1 + rgamma( n=S, shape=1, scale=0.75))
    }
  }

  if(family %in% c('bernoulli','binomial')) link <- make.link('logit')
  if(family %in% c('poisson','ippm','negative.binomial','tweedie')) link <- make.link('log')
  if(family %in% c('gaussian')) link <- make.link('identity')
  if(family %in% 'ippm') {
    grid <- simulate_ippm_grid(X,W)
    grid2D <- grid$grid2D
    X <- grid$X
    W <- grid$W
  }

  ## simulate the groups and fitted values.
  fitted <- matrix(0, dim(X)[1], S)
  group <- rep(0, S)
  if(npu>0) {
    eta_all <- U %*% delta
  } else {
    eta_all <- rep(0,nrow(X))
  }
  for (ss in seq_len(S)) {
    gg <- ceiling(stats::runif(1) * G)
    eta_spp <- W %*% c(alpha[ss],gamma[ss,])
    eta_mix <- X %*% beta[gg, ]
    eta <- eta_spp + eta_mix + eta_all + offset
    fitted[, ss] <- link$linkinv(eta)
    group[ss] <- gg
  }

  if( family=="bernoulli")
    outcomes <- matrix(rbinom(n * S, 1, as.numeric( fitted)), nrow = n, ncol = S)
  if( family=="binomial")
    outcomes <- matrix(rbinom(n * S, rep(size, S), as.numeric( fitted)), nrow = n, ncol = S)
  if( family=="poisson")
    outcomes <- matrix(rpois(n * S, lambda=as.numeric( fitted)), nrow = n, ncol = S)
  if( family=="ippm")
    outcomes <- simulate_ippm_outcomes(X, W, S, grid2D, fitted)
  if( family=="negative.binomial")
    outcomes <- matrix(rnbinom(n * S, mu=as.numeric( fitted), size=1/rep(exp(logTheta), each=n)), nrow = n, ncol = S)
  if (family=="tweedie")
    outcomes <- matrix(fishMod::rTweedie(n * S, mu = as.numeric(fitted),
                                         phi = rep(exp(logTheta), each = n),
                                         p = rep(powers, each = n)), nrow = n, ncol = S)
  if( family=="gaussian")
    outcomes <- matrix( rnorm( n=n*S, mean=as.numeric( fitted), sd=rep( exp(logTheta), each=n)), nrow=n, ncol=S)

  pi <- tapply(group, group, length)/S

  if (family=='ippm'){
    res <- outcomes$mm
    wts <- outcomes$weights
    colnames(wts) <- all.vars(sam_org)[1:S]
    colnames(fitted) <- all.vars(sam_org)[1:S]
  } else {
    colnames(fitted) <- all.vars(sam_org)[1:S]
    colnames(outcomes) <-  all.vars(sam_org)[1:S]
    if(ncol(W)>1){
      res <- data.frame(outcomes,const=1, X, W[,-1,drop=FALSE])
    } else {
      res <- data.frame(outcomes,const=1, X)
    }
    if(npu>0) res <- cbind(res,U)
  wts <- NULL
  }
  attr(res, "SAMs") <- group
  attr(res, "pis") <- pi
  attr(res, "alpha") <- alpha
  attr(res, "beta") <- beta
  attr(res, "gamma") <- gamma
  attr(res, "delta") <- delta
  attr(res, "logTheta") <- logTheta
  attr(res, "power") <- powers
  attr(res, "mu") <- fitted
  attr(res, "ippm_weights") <- wts
  attr(res, "size") <- size
  attr(res, "offset") <- offset
  return(res)
}


###### internal fitting functions #####
"fit.ecm.sam" <- function(y, X, W, U=NULL, spp_weights, site_spp_weights, offset,
                          y_is_na, G, S, disty, size, powers, control, starting.sam = NULL){

  n <- nrow(y)
  starting.sam <-  NULL
  bestOfAllMods <- list( logl=-Inf)
  for(t in seq_len(control$ecm_refit)) {
    if(!control$quiet) message(paste0("ECM restart ",t," of ", control$ecm_refit,""))
    starting.sam <- get_initial_values_sam(y = y, X = X, W = W, U = U,
                                           site_spp_weights = site_spp_weights,
                                           offset = offset, y_is_na = y_is_na,
                                           G = G, S = S,
                                           disty = disty, size = size,
                                           control = control)


    fits <- list()
    taus <- starting.sam$taus
    fits$pis <- colMeans(taus)
    fits$alpha <- starting.sam$alpha
    fits$betas <- starting.sam$beta
    fits$gamma <- starting.sam$gamma
    fits$delta <- starting.sam$delta
    fits$theta <- starting.sam$theta
    # if(fits$theta)

    diff.logl <- -100; cw.logl <- -Inf; new.logl <- 10
    mod <- list(logl = -Inf)

    counter <- 1
    while((diff.logl < 1-control$ecm_reltol | diff.logl > 1+control$ecm_reltol) & counter <= control$ecm_steps) {
      # if(diff.logl > 1) break;
      ## If any pi hits 0, restart with new random starts
    if(any(fits$pis < 0.0005)) {
    starting.sam <- get_initial_values_sam(y = y, X = X, W = W, U = U,
                                           site_spp_weights = site_spp_weights,
                                           offset = offset, y_is_na = y_is_na,
                                           G = G, S = S,
                                           disty = disty, size = size,
                                           control = set_control_sam(list(init_method = "random2")))

      fits <- list()
      taus <- starting.sam$taus
      fits$pis <- colMeans(taus)
      fits$alpha <- starting.sam$alpha
      fits$betas <- starting.sam$beta
      fits$gamma <- starting.sam$gamma
      fits$delta <- starting.sam$delta
      fits$theta <- starting.sam$theta
      }

      ## CM-step
      ## optimise the spp coefs - alpha & gamma if present
      fm_sppParam <- plapply(seq_len(S), apply_optimise_sppParam,
                             y, X, W, U, taus, fits,
                             site_spp_weights, offset, y_is_na,
                             disty, size, powers,
                             .parallel = control$cores,
                             .verbose = FALSE)
      new.sp.params <- do.call(rbind,lapply(fm_sppParam, `[[`, 1))
      fits$alpha <- update_coefs(fits$alpha,new.sp.params[,1])
      if(ncol(W)>1){
        fits$gamma <- update_coefs(fits$gamma,new.sp.params[,-1])
      } else {
        fits$gamma <- rep(-99999,S)
      }

      ## update the dispersion parameters.
      if(disty %in% c(4,5,6)) {
        fm_theta <- plapply(seq_len(S), apply_optimise_thetaParams,
                            y, X, W, U, taus, fits,
                            site_spp_weights, offset, y_is_na,
                            disty, size, powers,
                            .parallel = control$cores,
                            .verbose = FALSE)
        new.sp.thetas <- unlist(lapply(fm_theta, `[[`, 1))
        fits$theta <- update_coefs(fits$theta,new.sp.thetas)
      } else {
        fits$theta <- rep(-99999,S)
      }

      ## Update the all parameter.
      if(!is.null(U)){
        new.deltas <- stats::optimize(f = llogl.allParams, interval = c(-30,30), maximum = TRUE,
                               y = y, X = X, W = W, U=U, taus = taus, fits = fits,
                               site_spp_weights = site_spp_weights,
                               offset = offset, y_is_na = y_is_na,
                               disty = disty, size = size, powers = powers)$maximum
        fits$delta <- update_coefs(fits$delta, new.deltas)
      }

      ## Update archetype-specific coefficients
      fm_beta <- plapply(seq_len(G), apply_optimise_betas,
                         y, X, W, U, site_spp_weights, offset,
                         y_is_na, disty, taus, fits, size, powers,
                         .parallel = control$cores,
                         .verbose = FALSE)
      new.betas <- do.call(rbind,fm_beta)
      fits$beta <- update_coefs(fits$beta,new.betas)

      ## E-step
      do.estep <- try(e.step(y, X, W, U, site_spp_weights,
                             offset, y_is_na, disty,
                             fits, size, powers),
                      silent=TRUE)
      if(!is.finite(do.estep$logl)) {
        mod <- list(logl = -Inf); break;
      }
      if(is.finite(do.estep$logl)) {
        taus <- do.estep$taus
        new.logl <- do.estep$logl
        diff.logl <- abs(new.logl/cw.logl)

        if(!control$quiet) cat("Iteration: ", counter, "| New loglik", round(new.logl,3),  "| Ratio loglik", diff.logl, "\n");

        cw.logl <- new.logl
        counter <- counter + 1
      }
    }

    if(new.logl > mod$logl) {
      # mod$S <- S; mod$G <- G; mod$npx <- ncol(X); mod$npw <- ifelse(ncol(W)>1,ncol(W),0);
      # mod$npu <- ifelse(!is.null(U),ncol(U),0); mod$n <- n; mod$disty <- disty;
      mod <- list(logl = new.logl, alpha = fits$alpha, beta = fits$beta,
                  eta = additive_logistic(fits$pis, inv = TRUE)[-G],
                  gamma = fits$gamma, delta = fits$delta,
                  theta = fits$theta, pis = fits$pis, taus = taus)
      mod$diff.logl <- diff.logl
      mod$counter <- counter
      names(mod$pis) <- paste("Archetype", 1:G, sep = "");
      mod$taus <- as.data.frame(round(mod$taus,5)); colnames(mod$taus) <- names(mod$pis); rownames(mod$taus) <- paste("Sp",1:S,sep="")
      mod$entropy <- -sum(as.vector(mod$taus[mod$taus != 0])*log(unlist(mod$taus[mod$taus != 0])));
      mod$beta <- round(mod$beta,6); rownames(mod$beta) <- names(mod$pis);
      names(mod$alpha) <- rownames(mod$taus)
      if(ncol(W)>1) {
        rownames(gamma) <- rownames(mod$taus)
        colnames(gamma) <- colnames(W[,-1,drop=FALSE])
      } else {
        names(mod$gamma) <- rownames(mod$taus)
      }
      if(!is.null(U)) names(mod$beta) <- colnames(U)
      names(mod$theta) <- rownames(mod$taus)

      mod$effect.param <- length(mod$beta) + G + S + ifelse(ncol(W)>1,ncol(W)-1,0) + ifelse(!is.null(U),ncol(U),0) + ifelse(disty %in% c(4,5,6),S,0)
      mod$ics <- -2*mod$logl + c(2,log(S),log(S))*mod$effect.param + c(0,0,2*mod$entropy)
      names(mod$ics) <- c("AIC","BIC with log(# of species)", "ICL")
      mod$n <- n
    }
    if( mod$logl > bestOfAllMods$logl)
      bestOfAllMods <- mod
  }

  return(bestOfAllMods)
}

"e.step" <- function(y, X, W, U, site_spp_weights, offset,
                     y_is_na, disty,
                     fits, size, sp.powers,
                     get.fitted = FALSE) {
  S <- ncol(y); n <- nrow(y); G <- length(fits$pis)
  out.taus <- matrix(0,S,G)
  sp.logl <- rep(0,S) ## Species specific incomplete logL

  if(disty%in%c(1,7)) link <- make.link('logit')
  if(disty%in%c(2,3,4,5)) link <- make.link('log')
  if(disty%in%6) link <- make.link('identity')

  if(get.fitted) fitted.values <- array(0,dim=c(G,n,S))

  if(!is.null(U)) all.etas <- as.matrix(U)%*%c(fits$delta)
  else all.etas <- rep(0,n)

  for(gg in seq_len(G)) {
    mix.etas <- as.matrix(X)%*%fits$beta[gg,]

    for(ss in seq_len(S)) {
      sp_idx <- !y_is_na[,ss]
      if(ncol(W)>1){
        spp.etas <- as.matrix(W) %*% c(fits$alpha[ss],fits$gamma[ss,])
      }else{
        spp.etas <- as.matrix(W) %*% c(fits$alpha[ss])
      }
      new.etas <- all.etas + mix.etas + spp.etas + offset
      if(disty %in% 1) {
        check.p <- link$linkinv(new.etas)
        check.p[check.p < 1e-4] <- 1e-4; check.p[check.p > (1-1e-4)] <- (1-1e-4);
        if(get.fitted) fitted.values[gg,,ss] <- check.p
        out.taus[ss,gg] <- (sum(dbinom(y[,ss], 1, prob = check.p, log = TRUE)))
      }
      if(disty %in% 2) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.taus[ss,gg] <- (sum(dpois(y[,ss], lambda = link$linkinv(new.etas), log = TRUE)))
      }
      if(disty %in% 3) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.taus[ss,gg] <- (y[sp_idx,ss] %*% new.etas - site_spp_weights[sp_idx,ss] %*% link$linkinv(new.etas))
      }
      if(disty %in% 4) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.taus[ss,gg] <- (sum(dnbinom(y[,ss], mu = link$linkinv(new.etas), size = 1/fits$theta[ss], log = TRUE)))
      }
      if(disty %in% 5) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.taus[ss,gg] <- (sum(fishMod::dTweedie(y[,ss], mu = link$linkinv(new.etas), phi = fits$theta[ss], p = sp.powers[ss], LOG = TRUE)))
      }
      if(disty %in% 6) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.taus[ss,gg] <- (sum(dnorm(y[,ss], mean = link$linkinv(new.etas), sd = sqrt(fits$theta[ss]), log = TRUE)))
      }
      if(disty %in% 7) {
        check.p <- link$linkinv(new.etas)
        check.p[check.p < 1e-4] <- 1e-4; check.p[check.p > (1-1e-4)] <- (1-1e-4);
        if(get.fitted) fitted.values[gg,,ss] <- check.p
        out.taus[ss,gg] <- (sum(dbinom(y[,ss], size, prob = check.p, log = TRUE)))
      }

    }
  }

  for(ss in seq_len(S)) {
    eps <- max(out.taus[ss,])
    sp.logl[ss] <- log(sum(fits$pis*exp(out.taus[ss,]-eps))) + eps
  }
  for(ss in seq_len(S)) {
    for(gg in seq_len(G)) {
      out.taus[ss,gg] <- exp((log(fits$pis[gg]) + out.taus[ss,gg]) - sp.logl[ss])
    }
  }

  full.logl <- sum(sp.logl)
  out.list <- list(taus = out.taus, sp.logl = sp.logl, logl = full.logl)
  if(get.fitted) out.list$fitted = fitted.values
  return(out.list)
}

"apply_optimise_sppParam" <- function(ss, y, X, W, U, taus, fits,
                                      site_spp_weights, offset, y_is_na,
                                      disty, size, powers){
  if(disty%in%6) intercept.range <- c(-200,200)
  else intercept.range <- c(-30,30)
  update.sppParams <- stats::optimize(f = llogl.sppParams, interval = intercept.range, maximum = TRUE, ss = ss,
                               y = y, X = X, W = W, U=U, taus = taus, fits = fits,
                               site_spp_weights = site_spp_weights,
                               offset = offset, y_is_na = y_is_na,
                               disty = disty, size = size, powers = powers)$maximum
  return(update.sppParams)
}

"apply_optimise_thetaParams" <- function(ss, y, X, W, U, taus, fits,
                                   site_spp_weights, offy, y_is_na,
                                   disty, size, powers){
  update.theta <- suppressWarnings(stats::optimize(f = llogl.thetaParams, interval = c(0.001,100),
                                                   maximum = TRUE, ss = ss,
                           y = y, X = X, W = W, U=U, taus = taus, fits = fits,
                           site_spp_weights = site_spp_weights,
                           offset = offy, disty = disty, size = size, powers = powers)$maximum)
  return(update.theta)
}

"apply_optimise_betas" <- function(gg, y, X, W, U, site_spp_weights, offset,
                                   y_is_na, disty, taus, fits, size, powers){

  n <- nrow(y)
  if(disty %in% c(2,3,4,5)){
  glmnet.family <- "poisson"; glm.family <- poisson();
  }
  if(disty %in% c(1,7)){
  glmnet.family <- "binomial"; glm.family <- binomial();
  }
  if(disty %in% c(6)){
  glmnet.family <- "gaussian"; glm.family <- gaussian();
  }
        if(disty %in% 4) get.mus <- e.step(y, X, W, U, site_spp_weights,
                                           offset, y_is_na, disty,
                                           fits, size, powers,  get.fitted = TRUE)$fitted

        Y_s <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
        size_s <- matrix(rep(size,ncol(y)),nrow(y),ncol(y))[!y_is_na]
        X_no_NA <- list()
        for (jj in 1:ncol(y)){
          X_no_NA[[jj]] <- X[!y_is_na[,jj],,drop=FALSE]
        }
        X_s <- do.call(rbind, X_no_NA)

        if(!is.null(U)){
          U_no_NA <- list()
          for (jj in 1:ncol(y)){
            U_no_NA[[jj]] <- U[!y_is_na[,jj],,drop=FALSE]
          }
          U_s <- do.call(rbind, U_no_NA)
          offy3 <- as.matrix(U_s) %*% fits$delta
        } else {
          offy3 <- rep(0,nrow(X_s))
        }

        n_ys <- sapply(X_no_NA,nrow)
        tau.weights <- rep(taus[,gg,drop=FALSE],c(n_ys))
        if(disty %in% 4) tau.weights <- rep(taus[,gg,drop=FALSE],c(n_ys))/(1+rep(fits$theta,each=n)*as.vector(get.mus[gg,,]))
        site.weights <- as.matrix(as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
        obs.weights <- as.vector(tau.weights*site.weights)

        offy_mat <- replicate(ncol(y),offset)
        offy1 <- as.matrix(unlist(as.data.frame(offy_mat[!y_is_na])))
        if(ncol(W)>1){
          offy2 <- as.matrix(W) %*% t(cbind(fits$alpha,fits$gamma))
        } else {
          offy2 <- as.matrix(W) %*% c(fits$alpha)
        }
        offy2 <- as.matrix(unlist(as.data.frame(offy2[!y_is_na])))
        offy <- as.numeric(offy1 + offy2 + offy3)

        # which family to use?
        if (disty==3){
          Y_s <- as.vector(Y_s/site.weights)
        }
        if(disty%in%c(7)){
          Y_s <- as.matrix(cbind(Y_s,size_s-Y_s))
        }

        if(ncol(X_s)==1) X_s <- cbind(1,X_s)

        if (disty==7){
          fit1 <- try(glm2::glm.fit2(y=Y_s, x=as.data.frame(cbind(1,X_s)),weights=c(obs.weights+1e-6),
                               family=glm.family, offset=offy), silent=FALSE)
          if (any(class(fit1)[1] %in% 'try-error')){
            new.betas <- rep(NA, ncol(X))
            names(new.betas) <- colnames(cbind(X,W))
          } else {
             new.betas <- coef(fit1)[-1]
          }
        } else {
          fit1 <- glmnet::glmnet(x = as.matrix(X_s), y = Y_s, family = glmnet.family, weights = obs.weights+1e-6, offset = offy, nlambda = 100, intercept = FALSE)
        new.betas <- coef(fit1)[,ncol(coef(fit1))][-1]
        }
        # print(new.betas)
        if(ncol(X)==1)new.betas <- new.betas[-1]

        return(new.betas)

}

"llogl.sppParams" <- function(x, ss, y, X, W, U, taus, fits, site_spp_weights, offset, y_is_na, disty, size, powers) {

  S <- ncol(y); n <- nrow(y); G <- ncol(taus)
  out <- 0
  sp_idx <- !y_is_na[,ss]

  if(disty%in%c(1,7)) link <- make.link('logit')
  if(disty%in%c(2,3,4,5)) link <- make.link('log')
  # if(disty%in%6) link <- make.link('identity')

  if(!is.null(U)) all.etas <- as.matrix(U)%*%c(fits$delta)
  else all.etas <- rep(0,n)

  mix.etas <- as.matrix(X)%*%t(fits$betas)

  for(gg in seq_len(G)) {
    eta <- as.matrix(W)%*%x + mix.etas[,gg] + all.etas + offset
    if(disty %in%  1)
      out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = 1, prob =  link$linkinv(eta), log = TRUE))
    if(disty %in%  2)
      out <- out + sum(taus[ss,gg]*dpois(y[,ss], lambda =  link$linkinv(eta), log = TRUE))
    if(disty %in%  3)
      out <- out + sum(taus[ss,gg]*(y[sp_idx,ss] %*% eta - site_spp_weights[sp_idx,ss] %*%  link$linkinv(eta)))
    if(disty %in%  4)
      out <- out + sum(taus[ss,gg]*dnbinom(y[,ss], mu =  link$linkinv(eta), size = 1/fits$theta[ss], log = TRUE))
    if(disty %in%  5)
      out <- out + sum(taus[ss,gg]*fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
    if(disty %in%  6)
      out <- out + sum(taus[ss,gg]*dnorm(y[,ss], mean =  eta, sd = sqrt(fits$theta[ss]), log = TRUE))
    if(disty %in%  7)
      out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = size, prob = link$linkinv(eta), log = TRUE))
  }
  return(out)
}


"llogl.thetaParams" <- function(x, ss, y, X, W, U, taus, fits, site_spp_weights, offset, disty, size, powers) {

  S <- ncol(y); n <- nrow(y); G <- ncol(taus)
  out <- 0
  if(disty%in%c(4,5)) link <- make.link('log')
  if(disty%in%6) link <- make.link('identity')


  if(!is.null(U)) all.etas <- as.matrix(U)%*%c(fits$delta)
  else all.etas <- rep(0,n)

  if(ncol(W)>1){
    spp.etas <- as.matrix(W) %*% c(fits$alpha[ss],fits$gamma[ss,])
  }else{
    spp.etas <- as.matrix(W) %*% c(fits$alpha[ss])
  }
  mix.etas <- as.matrix(X) %*% t(fits$beta)
  for(gg in seq_len(G)){
    eta <- spp.etas + mix.etas[,gg] + all.etas + offset
   if(disty==4)
     out <- out + sum(taus[ss,gg]*dnbinom(y[,ss], mu = link$linkinv(eta), size = 1/x, log = TRUE));
   if(disty==5)
     out <- out + sum(taus[ss,gg]*fishMod::dTweedie(y[,ss], mu = link$linkinv(eta), phi = x, prob = powers[ss], LOG = TRUE));
   if(disty==6)
     out <- out + sum(taus[ss,gg]*dnorm(y[,ss], mean = eta, sd = sqrt(x), log = TRUE))
  }
  return(out)
}



"llogl.allParams" <- function(x, y, X, W, U, taus, fits, site_spp_weights, offset, y_is_na, disty, size, powers) {

  S <- ncol(y); n <- nrow(y); G <- ncol(taus)
  out <- 0

  if(disty%in%c(1,7)) link <- make.link('logit')
  if(disty%in%c(2,3,4,5)) link <- make.link('log')
  if(disty%in%6) link <- make.link('identity')

  all.etas <- as.matrix(U)%*%x

  mix.etas <- as.matrix(X)%*%t(fits$beta)

  for (ss in seq_len(S)){

    sp_idx <- !y_is_na[,ss]
    if(ncol(W)>1){
      spp.etas <- as.matrix(W) %*% c(fits$alpha[ss],fits$gamma[ss,])
    }else{
      spp.etas <- as.matrix(W) %*% c(fits$alpha[ss])
    }

    for(gg in seq_len(G)) {
      eta <- spp.etas + mix.etas[,gg] + all.etas + offset
      if(disty %in%  1)
        out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = 1, prob =  link$linkinv(eta), log = TRUE))
      if(disty %in%  2)
        out <- out + sum(taus[ss,gg]*dpois(y[,ss], lambda =  link$linkinv(eta), log = TRUE))
      if(disty %in%  3)
        out <- out + sum(taus[ss,gg]*(y[sp_idx,ss] %*% eta - site_spp_weights[sp_idx,ss] %*%  link$linkinv(eta)))
      if(disty %in%  4)
        out <- out + sum(taus[ss,gg]*dnbinom(y[,ss], mu =  link$linkinv(eta), size = 1/fits$theta[ss], log = TRUE))
      if(disty %in%  5)
        out <- out + sum(taus[ss,gg]*fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], prob = powers[ss], LOG = TRUE))
      if(disty %in%  6)
        out <- out + sum(taus[ss,gg]*dnorm(y[,ss], mean =  (eta), sd = sqrt(fits$theta[ss]), log = TRUE))
      if(disty %in%  7)
        out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = size, prob = link$linkinv(eta), log = TRUE))
    }

  }
  return(out)
}


"starting_values_wrapper" <- function(y, X, W, U, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, size, powers, control){
  if(control$ecm_prefit){
    if(!control$quiet)message('Using ECM algorithm to find starting values; using ',
                              control$ecm_refit,'refits')
    emfits <- fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights,
                            offset, y_is_na, G, S, disty, size, powers, control)
    start_vals <- list(alpha = (emfits$alpha),
                       beta = (emfits$beta),
                       gamma = (emfits$gamma),
                       delta = (emfits$delta),
                       theta = (emfits$theta),
                       pis = (emfits$pis))
  } else {
    message('If you are choosing to not use the ECM algorithm to estimate starting values\nwe recommend that you at least run a multifits to optimise the loglikelihood, see "species_mix.multifit".')
    if(!control$quiet)message('You are not using the EM algorith to find
starting values;\n starting values are generated using ',control$init_method,
                              '.')
    starting_values <- get_initial_values_sam(y = y, X = X, W= W, U = U,
                                              site_spp_weights = site_spp_weights,
                                              offset = offset, y_is_na = y_is_na,
                                              G = G, S = S,
                                              disty=disty,size=size, powers = powers,
                                              control = control)
    start_vals <- list(alpha=(starting_values$alpha),
                       beta=(starting_values$beta),
                       gamma=(starting_values$gamma),
                       delta=(starting_values$delta),
                       theta=(starting_values$theta),
                       pis=(starting_values$pis))
  }

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pis, inv = TRUE)[-G]
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)
  #put dispersion parameters on the log-scale for c++
  if(disty%in%c(4,5,6)) start_vals$theta <- transform.theta(start_vals$theta, disty=disty, max.theta = control$max.theta)
  return(start_vals)
}

## function for starting values using penalities
"get_initial_values_sam" <- function(y, X, W, U=NULL, site_spp_weights,
                                   offset, y_is_na, G, S, disty, size, powers, control) {

  n <- nrow(y)
  prev.min <- floor(n*control$minimum_sites_occurrence);
  sel.omit.spp <- which(colSums(y>0) <= prev.min)
  if(length(sel.omit.spp)==0) sel.omit.spp <- -1*(1:S)

  starting.sam <- list(alpha = rep(0,S), theta = rep(1,S));
  if(!control$quiet) message("Initialising starting values")


  if(is.null(U)){
    datX <- cbind(W,X)
  } else {
    datX <- cbind(W,X,U)
  }
  tmpform <- as.formula(paste0("mvabund::mvabund(y) ~ - 1 + ",paste0(colnames(datX),collapse = "+")))

  if(disty%in%1){ # bernoulli
    fit1 <- mvabund::manyglm(tmpform, data = datX, offset = offset, family = "binomial")
    coefs <- t(fit1$coefficients)
    starting.sam$theta <- rep(-99999,S)
  }
  if(disty%in%2){ # poisson
    fit1 <- mvabund::manyglm(tmpform, data = datX, offset = offset, family = "poisson")
    coefs <- t(fit1$coefficients)
    starting.sam$theta <- rep(-99999,S)
  }
  if(disty%in%3){ # ippm
    fit1 <- many.fit(y, X, W, U, site_spp_weights,
                      offset, y_is_na, G, S, disty, size, powers, control)
    coefs <- fit1$coefficients
  }
  if(disty%in%4){ # negative binomial
    fit1 <- mvabund::manyglm(tmpform, data = datX,offset = offset, family = "negative.binomial")
    coefs <- t(fit1$coefficients)
    starting.sam$theta <- ifelse(fit1$phi<1,1,fit1$phi)
  }
  if(disty%in%5){ # tweedie
    fit1 <- many.fit(y, X, W, U, site_spp_weights,
                     offset, y_is_na, G, S, disty,
                     size, powers, control)
    starting.sam$theta <- exp(fit1$theta)
    coefs <- fit1$coefficients
  }
  if(disty%in%6){ # gaussian
    fit1 <- suppressWarnings(mvabund::manylm(tmpform,data = datX))
    starting.sam$theta <- (stats::deviance(fit1)/fit1$df.residual)
    # log( sqrt( colSums((y - fit1$fitted.values)^2)/nrow(y)))
    coefs <- t(fit1$coefficients)
  }
  if(disty%in%7){  # binomial - with variable size
    fit1 <- many.fit(y, X, W, U, site_spp_weights,
                     offset, y_is_na, G, S, disty, size, powers, control)
    coefs <- fit1$coefficients
    starting.sam$theta <- rep(-99999,S)
  }

  # a bit of book keeping
  covarNames <- colnames(datX)
  alphaIdx <- match(colnames(W[,1,drop=FALSE]),covarNames)
  betaIdx <- match(colnames(X),covarNames)

  # alphas/intercepts
  starting.sam$alpha <- coefs[,alphaIdx, drop=FALSE]

  # betas to setup archetype parameters
  spp.beta <- coefs[-sel.omit.spp,betaIdx, drop=FALSE]

  # gammas for species intercepts
  if(ncol(W)>1){
    gammaIdx <- match(colnames(W[,-1,drop=FALSE]),covarNames)
    starting.sam$gamma <- coefs[,gammaIdx, drop=FALSE]
  } else {
    starting.sam$gamma <- -99999
  }

  # delta for entire dataset
  if(!is.null(U)){
    deltaIdx <- match(colnames(U),covarNames)
    starting.sam$delta <- colMeans(coefs[-sel.omit.spp,deltaIdx,drop=FALSE])
  } else {
    starting.sam$delta <- -99999
  }


  if(G==1) control$init_method <- 'kmeans'

  if(control$init_method=='kmeans'){
    if(!control$quiet) message("Initial groups parameter estimates by K-means clustering")
    fmmvnorm <- stats::kmeans(spp.beta, centers=G, iter.max = 200, nstart = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(spp.beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)
    colnames(grp_coefs) <- covarNames[betaIdx]
    rownames(grp_coefs) <-  paste("Archetype", 1:G, sep = "");
  }

  if(control$init_method=='kmed'){
    if(!control$quiet) message("Initial groups parameter estimates by K-medoids")
    mrwdist <- kmed::distNumeric(spp.beta, spp.beta, method = "mrw")
    fmmvnorm <- kmed::fastkmed(mrwdist, ncluster = G, iterate = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- spp.beta[fmmvnorm$medoid,,drop=FALSE]
    grp_coefs <- matrix(grp_coefs,nrow=G)
    colnames(grp_coefs) <- covarNames[betaIdx]
    rownames(grp_coefs) <-  paste("Archetype", 1:G, sep = "");
  }

  if(control$init_method=='random2'){
    if(!control$quiet) message("Initial groups parameter estimates by K-means clustering with random noise")
    fmmvnorm <- stats::kmeans(spp.beta, centers=G, nstart=50, iter.max = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(spp.beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G) # not check this...
    colnames(grp_coefs) <- covarNames[betaIdx]
    rownames(grp_coefs) <-  paste("Archetype", 1:G, sep = "");
    random_coefs <- sam_random_inits(alpha = starting.sam$alpha,beta = grp_coefs,
                                              gamma = starting.sam$gamma, delta = starting.sam$delta,
                                              theta = starting.sam$theta, S, G, X, W, U,
                                              disty, mult=0.3, control$init_sd)
    starting.sam$alpha <- random_coefs[[1]]
    grp_coefs <- random_coefs[[2]]
    starting.sam$gamma <- random_coefs[[3]]
    starting.sam$delta <- random_coefs[[4]]
    if(disty%in%c(4,5,6)) starting.sam$theta <- random_coefs[[5]]
  }

  #get taus as starting values
  if(G==1){
    taus <- matrix(1,nrow=ncol(y), ncol = G)
  } else {
    taus <- matrix(0,nrow=ncol(y), ncol = G)
    if(length(sel.omit.spp)>0){
      for(j in 1:length((1:S)[-sel.omit.spp]))
        taus[(1:S)[-sel.omit.spp][j],fmmvnorm$cluster[j]] <- 1
      taus[sel.omit.spp,] <- matrix(runif(length(sel.omit.spp)*G),length(sel.omit.spp), G)
    } else {
      for(j in seq_len(S))taus[j,fmmvnorm$cluster[j]] <- 1
    }
  }

  taus <- taus/rowSums(taus)
  # starting.sam$alpha
  starting.sam$beta <- grp_coefs
  starting.sam$taus <- shrink_taus(taus, G=G)
  starting.sam$pis <- colMeans(taus)

  return(starting.sam)
}

"many.fit" <- function(y, X, W, U, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control){

  options(warn = -1)
  fm_sp_mods <-  plapply(seq_len(S), apply_species_fits, y, X, W, U,
                                  site_spp_weights, offset, y_is_na, disty, size, powers,
                                  .parallel = control$cores, .verbose = FALSE)

  alpha <- unlist(lapply(fm_sp_mods, `[[`, 1))
  alphaName <- "Intercept"
  if(ncol(X)==1){
    beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))[,-1,drop=FALSE]
  } else {
    beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))
  }
  betaNames <- colnames(X)
  if(ncol(X)==0){
    beta <- rep(-999999,G)
    betaNames <- "betaNAN!"
  }
  if(ncol(W)>1){
    gamma <- do.call(rbind,lapply(fm_sp_mods, `[[`, 3))
    gammaNames <- colnames(W[-1,,drop=FALSE])
  } else {
    gamma <- unlist(lapply(fm_sp_mods, `[[`, 3))
    gammaNames <- "gammaNAN!"
  }

  if(!is.null(U)){
    delta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 4))
    deltaNames <- colnames(U)
  } else {
    delta <- unlist(lapply(fm_sp_mods, `[[`, 4))
    deltaNames <- "deltaNAN!"
  }

  theta <- unlist(lapply(fm_sp_mods, `[[`, 5))

  res <- list()
  coefficients <- cbind(alpha,beta,gamma,delta)
  colnames(coefficients) <- c(alphaName,betaNames,gammaNames,deltaNames)
  rownames(coefficients) <- colnames(y)

  dropCovar <- grep("NAN!",colnames(coefficients))

  res$coefficients <- coefficients[,-dropCovar,drop=FALSE]
  res$theta <- theta
  return(res)
}


"apply_species_fits" <- function(ss, y, X, W, U = NULL, site_spp_weights,
                                     offset, y_is_na, disty, size, power){

  # which family to use?
  if(disty %in% c(1,7)) #binomials
    fam <- "binomial" #glmnet
  if(disty %in% c(2,3,4))
    fam <- "poisson"
  if(disty %in% 6)
    fam <- "gaussian"

  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.matrix(y[ids_i,ss])
  }

  if(ncol(X)==1){
    X<- cbind(1,X[ids_i,,drop=FALSE])
  }

  if(ncol(W) > 1){
    df <- cbind(X[ids_i,,drop=FALSE],W[ids_i,-1,drop=FALSE])
  } else {
    df <-   X[ids_i,,drop=FALSE]
  }

  if(!is.null(U)) df <- cbind(df,U[ids_i,,drop=FALSE])

  if(!disty %in% c(5,7)){

    lambda.seq <- sort(unique( c( seq( from=1/0.001, to=1, length=25),
                                   seq( from=1, to=.1, length=10))),
                       decreasing=TRUE)
    # if(disty==7) lambda.seq <- 0
    ft_sp <- try(glmnet::glmnet(y=outcomes, x=as.matrix(df),
                                family=fam, offset=offset[ids_i],
                                weights=as.numeric(site_spp_weights[ids_i,ss]),
                                alpha=0,
                                lambda=lambda.seq,
                                standardize=FALSE,
                                intercept=TRUE), silent=FALSE)
    locat.s <- 1/1
    my.coefs <- glmnet::coef.glmnet(ft_sp, s=locat.s)
    if( any( is.na( my.coefs))){  #just in case the model is so badly posed that mild penalisation doesn't work...
      my.coefs <- glmnet::coef.glmnet(ft_sp, s=lambda.seq)
      lastID <- apply( my.coefs, 2, function(x) !any( is.na( x)))
      lastID <- tail( (seq_along( lastID))[lastID], 1)
      my.coefs <- my.coefs[,lastID]
    }
    if (any(class(ft_sp) %in% 'try-error')){
      my_coefs <- rep(NA, ncol(X[ids_i,]))
      names(my_coefs) <- colnames(cbind(X[ids_i,,drop=FALSE],W[ids_i,,drop=FALSE],U[ids_i,,drop=FALSE]))
    } else {
      if(ncol(X)==1) my_coefs <- t(as.matrix(my.coefs[-1]))
      my_coefs <- t(as.matrix(my.coefs))
    }
    theta <- -99999
    if( disty == 4){
      tmp <- MASS::theta.mm(outcomes, as.numeric(predict(ft_sp, s=locat.s,
                                                         type="response",
                                                         newx=as.matrix(df),
                                                         newoffset=offset[ids_i])),
                            weights=as.matrix(site_spp_weights[ids_i,ss]),
                            dfr=length(outcomes), eps=1e-4)
      if(tmp>5) tmp <- 5
      theta <- tmp
    }
    if( disty == 6){
      preds <- as.numeric( predict(ft_sp, s=locat.s, type="link",
                                   newx=as.matrix(df), newoffset=offset[ids_i]))
      #should be something like the resid standard
      theta <- log( sqrt( sum((outcomes - preds)^2)/length(outcomes)))
    }
  }

  if (disty==5) { #Tweedie needs an unconstrained fit.  May cause problems in some cases, especially if there is quasi-separation...
      df3 <- data.frame(y=outcomes, offy=offset, Intercept= 1, df)
      tmp.fm1 <- fishMod::tglm(y~-1+.-offy+offset( offy),
                               wts=c(site_spp_weights[,ss]),
                               data=df3, p=power[ss], vcov=FALSE,
                               residuals=FALSE, trace=0)
      my_coefs <- t(as.matrix(tmp.fm1$coef,ncol=1))
      theta <- log(tmp.fm1$coef["phi"])
  }

  if (disty==7){
      outcomes <- as.matrix(cbind(y[ids_i,ss],size[ids_i]-y[ids_i,ss]))
      fit_init <- glm2::glm.fit2(y=outcomes, x=cbind(1,df), family = binomial())
      my_coefs <- t(as.matrix(fit_init$coefficients))
      theta <- -99999
    }

  # species intercpets
  alpha <- my_coefs[1]
  # mixture coefs
  beta <- my_coefs[match(colnames(X), colnames(my_coefs))]
  # species coefs apart from intercept
  if(ncol(W)>1) gamma <-  my_coefs[match(colnames(W[,-1,drop=FALSE]), colnames(my_coefs))] else gamma <- -99999
  if(!is.null(U)) delta <-  my_coefs[match(colnames(U), colnames(my_coefs))] else delta <- -99999

  return(list(alpha = alpha, beta = beta, gamma = gamma, delta = delta, theta = theta))
}


"sam_optimise" <- function(y, X, W, U, offset, spp_weights, site_spp_weights, y_is_na,
                           S, G, disty, size, powers, start_vals, control){

  inits <- unname(c(start_vals$alpha, start_vals$beta, start_vals$eta, start_vals$gamma,
             start_vals$delta, start_vals$theta))
  npx <- as.integer(ncol(X))
  n <- as.integer(nrow(X))

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

  # sort of the species formula structures for cpp
  if(ncol(W)>1){
    npw <- as.integer(ncol(W[,-1,drop=FALSE]))
    control$optiPart <- as.integer(1)
    gamma.score <- as.numeric(matrix(NA, nrow=S, ncol=ncol(W)))
    Wcpp <- W[,-1,drop=FALSE]
  } else {
    npw <- as.integer(1)
    control$optiPart <- as.integer(0)
    gamma.score <- -99999
    Wcpp <- matrix(1,nrow = n, ncol=1)
  }

if(!is.null(U)) {
    Ucpp <- U
    npu <- as.integer(ncol(U))
    control$optiAll <- as.integer(1)
    delta.score <- as.numeric(matrix(NA, ncol=ncol(U)))
  } else {
    delta.score <- -99999
    control$optiAll <- as.integer(0)
    Ucpp <- matrix(1,nrow = n,ncol=1)
    npu <- as.integer(1) # a dummy variable to stop c++ issues.
  }

  if(disty%in%c(4,5,6)){
    control$optiDisp <- as.integer(1)
    theta.score <- as.numeric(rep(NA, length(theta)))
  }else{
    control$optiDisp <- as.integer(0)
    theta.score <- -99999
  }
  scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,
                                       gamma.score,delta.score,theta.score))))
  control$conv <- as.integer(0)

  #model quantities
  pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, S))
  loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

  if(control$print_cpp_start_vals)print_starting_values(as.integer(S),
                                                        as.integer(G),
                                                        as.integer(npx),
                                                        as.integer(npw),
                                                        as.integer(npu),
                                                        as.integer(n),
                                                        as.double(alpha),
                                                        as.double(beta),
                                                        as.double(gamma),
                                                        as.double(eta),
                                                        as.double(delta),
                                                        as.double(theta))

  #c++ call to optimise the model (needs pretty good starting values)
  tmp <- .Call("species_mix_cpp",
               as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
               as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
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
               theta.score, as.integer(control$getscores_cpp), as.numeric(scores),
               # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
               pis_out, mus, loglikeS, loglikeSG,
               # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
               as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
               as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
               as.integer(control$printparams_cpp),
               # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
               as.integer(control$optimise_cpp), as.integer(control$loglOnly_cpp),
               as.integer(control$derivOnly_cpp),
               # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
               PACKAGE = "ecomix")

  ret <- tmp
  ret$logl <- ret$logl * -1
  ret$mus <- array(mus, dim=c(n, S, G))

  ret$coefs <- list(alpha = ret$alpha,
                    beta = matrix(ret$beta,G,npx),
                    eta = ret$eta,
                    gamma = matrix(ret$gamma,S,npw),
                    delta = ret$delta,
                    theta = ret$theta)

  ret$names <- list(spp=colnames(y), SAMs=paste("Archetype.", seq_len(G), sep=""),
                    Xvars=colnames(X), Wvars=colnames(Wcpp), Uvars = colnames(Ucpp))

  ret$scores <- list(alpha.scores = alpha.score,
                     beta.scores = beta.score,
                     eta.scores=eta.score,
                     gamma.scores = gamma.score,
                     delta.scores=delta.score,
                     theta.scores=theta.score)

  ret$S <- S; ret$G <- G; ret$npx <- npx; ret$npw <- ifelse(ncol(W)>1,ncol(W),0);
  ret$npu <- ifelse(!is.null(U),ncol(U),0); ret$n <- n; ret$disty <- disty;
  ret$start.vals <- inits
  ret$loglikeSG <- matrix(loglikeSG,  nrow = S, ncol = G)  #for residuals
  ret$loglikeS <- loglikeS  #for residuals
  gc()
  return(ret)
}

"sam_optimise_tweedie" <- function(y, X, W, U, offset, spp_weights, site_spp_weights, y_is_na,
                                   S, G, disty, size, powers, start_vals, control){

    Tw.phi.func <- function( phi1, spp3){
      disp3 <- theta
      disp3[spp3] <- phi1
      tmp1 <- .Call("species_mix_cpp",
                   as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
                   as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
                   as.numeric(offset), as.numeric(spp_weights),
                   as.numeric(as.matrix(site_spp_weights)),
                   as.integer(as.matrix(!y_is_na)),
                   as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
                   as.integer(npw), as.integer(npu), as.integer(n),
                   as.integer(disty),as.integer(TRUE),
                   as.integer(control$optiPart),as.integer(control$optiAll),
                   # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
                   as.double(alpha), as.double(beta), as.double(eta),
                   as.double(gamma), as.double(delta), as.double(disp3),
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
                   theta.score, as.integer(control$getscores_cpp), as.numeric(scores),
                   # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
                   pis_out, mus, loglikeS, loglikeSG,
                   # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
                   as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
                   as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
                   as.integer(control$printparams_cpp),
                   # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
                   as.integer(FALSE), as.integer(TRUE), as.integer(FALSE),
                   # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
                   PACKAGE = "ecomix")
      return( -as.numeric( tmp1))
    }

    Tw.phi.func.grad <- function( phi1, spp3){
      disp3 <- theta
      disp3[spp3] <- phi1
      tmp.disp.score <- rep( -99999, S)

      tmp1 <- .Call("species_mix_cpp",
                    as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
                    as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
                    as.numeric(offset), as.numeric(spp_weights),
                    as.numeric(as.matrix(site_spp_weights)),
                    as.integer(as.matrix(!y_is_na)),
                    as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
                    as.integer(npw), as.integer(npu), as.integer(n),
                    as.integer(disty),as.integer(TRUE),
                    as.integer(control$optiPart),as.integer(control$optiAll),
                    # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
                    as.double(alpha), as.double(beta), as.double(eta),
                    as.double(gamma), as.double(delta), as.double(disp3),
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
                    theta.score, as.integer(control$getscores_cpp), as.numeric(scores),
                    # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
                    pis_out, mus, loglikeS, loglikeSG,
                    # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
                    as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
                    as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
                    as.integer(control$printparams_cpp),
                    # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
                    as.integer(FALSE), as.integer(FALSE), as.integer(TRUE),
                    # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
                    PACKAGE = "ecomix")

      return( -as.numeric( tmp.disp.score[spp3]))
    }

    inits <- unname(c(start_vals$alpha, start_vals$beta, start_vals$eta, start_vals$gamma,
                      start_vals$delta, start_vals$theta))
    npx <- as.integer(ncol(X))
    n <- as.integer(nrow(X))

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

    # sort of the species formula structures for cpp
    if(ncol(W)>1){
      npw <- as.integer(ncol(W[,-1,drop=FALSE]))
      control$optiPart <- as.integer(1)
      gamma.score <- as.numeric(matrix(NA, nrow=S, ncol=ncol(W)))
      Wcpp <- W[,-1,drop=FALSE]
    } else {
      npw <- as.integer(1)
      control$optiPart <- as.integer(0)
      gamma.score <- -99999
      Wcpp <- matrix(1,nrow = n, ncol=1)
    }

    if(!is.null(U)) {
      Ucpp <- U
      npu <- as.integer(ncol(U))
      control$optiAll <- as.integer(1)
      delta.score <- as.numeric(matrix(NA, ncol=ncol(U)))
    } else {
      delta.score <- -99999
      control$optiAll <- as.integer(0)
      Ucpp <- matrix(1,nrow = n,ncol=1)
      npu <- as.integer(1) # a dummy variable to stop c++ issues.
    }

    if(disty%in%c(4,5,6)){
      control$optiDisp <- as.integer(1)
      theta.score <- as.numeric(rep(NA, length(theta)))
    }else{
      control$optiDisp <- as.integer(0)
      theta.score <- -99999
    }
    scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,
                                         gamma.score,delta.score,theta.score))))
    control$conv <- as.integer(0)

    #model quantities
    pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
    mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
    loglikeS <- as.numeric(rep(NA, S))
    loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

    optimiseDisp <- FALSE
    kount <- 1
    tmp.new <- tmp.old <- -999999
    if( control$optimise_cpp){
      while( (abs( abs( tmp.new - tmp.old) / ( abs( tmp.old) + control$reltol)) > control$reltol | kount==1) & (kount < 15)){
        kount <- kount + 1
        tmp.old <- tmp.new
        message( "Updating Location Parameters: ", appendLF=FALSE)
        tmp <- .Call("species_mix_cpp",
                      as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
                      as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
                      as.numeric(offset), as.numeric(spp_weights),
                      as.numeric(as.matrix(site_spp_weights)),
                      as.integer(as.matrix(!y_is_na)),
                      as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
                      as.integer(npw), as.integer(npu), as.integer(n),
                      as.integer(disty),as.integer(TRUE),
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
                      theta.score, as.integer(control$getscores_cpp), as.numeric(scores),
                      # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
                      pis_out, mus, loglikeS, loglikeSG,
                      # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
                      as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
                      as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
                      as.integer(control$printparams_cpp),
                      # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
                      as.integer(control$optimise_cpp), as.integer(TRUE), as.integer(FALSE),
                      # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
                      PACKAGE = "ecomix")
        # tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
                     # as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
                     # alpha, tau, beta, gamma, disp, power,
                     # as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
                     # alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                     # pis, mus, logCondDens, logls,
                     # as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
                     # as.integer(control$optimise_cpp), as.integer(TRUE), as.integer( FALSE), as.integer(optimiseDisp), as.integer( FALSE), PACKAGE = "ecomix")
        message( "Updating Dispersion Parameters: ", appendLF=FALSE)
        for( ii in 1:S){
          tmp1 <- nlminb( theta[ii], Tw.phi.func, Tw.phi.func.grad, spp3=ii, control=list( trace=0))
          theta[ii] <- tmp1$par
          message( tmp1$objective, " ")
        }
        message( "")
        tmp.new <- -tmp1$objective
      }
    }
    tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
                 as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
                 alpha, tau, beta, gamma, disp, power,
                 as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
                 alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
                 pis, mus, logCondDens, logls,
                 as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
                 as.integer(FALSE), as.integer( TRUE), as.integer(TRUE), as.integer(TRUE), as.integer( FALSE), PACKAGE = "ecomix")

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

"transform.theta" <- function(thetas, disty, max.theta = NULL){

  if(disty%in%4){
    if(!is.null(max.theta)) {
      thetas <- ifelse(thetas>max.theta,max.theta,thetas)
    }
  }
  if(disty%in%5){

  }
  if(disty%in%6){
    thetas <- sqrt(thetas)
  }

  log.new.thetas <- log(thetas)

  return(log.new.thetas)
}


###### SAM internal functions ######

"additive_logistic" <- function (x,inv=FALSE) {
  if(inv){
    x <- log(x/x[length(x)])
    return(x)
  } else {

    x.t <- exp(x)
    x.t <- x.t/(1+sum(x.t))
    x.t[length(x.t)+1] <- 1-sum(x.t)
    return(x.t)
  }
}


"calc_info_crit_sam" <-  function(tmp) {
  k <- length(tmp$beta) + tmp$G + tmp$S + tmp$npw*tmp$S +
    tmp$npu + tmp$disty
  tmp$BIC <- -2 * tmp$logl + log(tmp$n) * k
  tmp$AIC <- -2 * tmp$logl + 2 * k
  return( tmp)
}

"calc_post_probs_sam" <- function( pis, logCondDens)  {
  logPostProbs <- log( pis) + logCondDens
  mset <- apply( logPostProbs, 1, max)
  logSums <- mset + log( rowSums( exp( logPostProbs-mset)))
  logPostProbs <- logPostProbs - logSums
  postProbs <- exp( logPostProbs)
  return( postProbs)
}

"check_species_formula" <- function(f){

  if(is.null(f))
    return(0)

  if(all(attr(stats::terms(f), 'intercept') == 1 & length(attr(stats::terms(f), 'factors')) == 0))
    return(1)

  if(all(attr(stats::terms(f), 'intercept') == 1 & length(attr(stats::terms(f), 'factors')) > 0))
    return(2)

}

"check_size_binomial" <- function(size, nsites) {
  if(is.null(size)) size <- rep(1,nsites)
  if(length(size)!=nsites)stop("Vector sites does not match then number of sites, if you are supplying these values please double check.")
  return(size)
}

"check_reponse_sam" <-function(outs) {
  nam <- colnames( outs)
  if( length( nam) == length( unique( nam)))
    return( length( nam))
  else
    return( FALSE)
}


"check_spp_weights" <- function(bb_weights, nS){

  if(!is.null(bb_weights)){
    if(all.equal(length(bb_weights),nS)) {
      spp_weights <- bb_weights
    } else {
      message('species weights for Bayesian bootstrap do not match the number of species in the model. Please check.')
      spp_weights <- rep(1,nS)
    }
  } else {
    spp_weights <- rep(1,nS)
  }
  return(spp_weights)
}

"clean_data_sam" <- function(data, form1, form2, form3=NULL, family){
  if(family=='ippm') na_rule <- "na.pass"
  else  na_rule <- "na.exclude"
  mf.X <- stats::model.frame(form1, data = data, na.action = na_rule)
  if(is.null(form3)){
    mf.U <- NULL
    if(!is.null( form2)){
      # deparse(form2)
      mf.W <- stats::model.frame(form2, data = data, na.action = na_rule)
      ids <- c( rownames( mf.W), rownames( mf.X))[duplicated( c( rownames( mf.W), rownames( mf.X)))]  #those rows of data that are good for both parts of the model.
      mf.X <- mf.X[rownames( mf.X) %in% ids,, drop=FALSE]
      mf.W <- mf.W[rownames( mf.W) %in% ids,, drop=FALSE]
    } else{
      mf.W <- NULL
      ids <- rownames( mf.X)
    }
  } else {
    if(!is.null( form2)){
      mf.W <- stats::model.frame(form2, data = data, na.action = na_rule)
      mf.U <- stats::model.frame(form3, data = data, na.action = na_rule)
      ids <- c( rownames( mf.W), rownames( mf.X), rownames( mf.U))[duplicated( c( rownames( mf.W), rownames( mf.X), rownames( mf.U)))]  #those rows of data that are good for both parts of the model.
      mf.X <- mf.X[rownames( mf.X) %in% ids,, drop=FALSE]
      mf.W <- mf.W[rownames( mf.W) %in% ids,, drop=FALSE]
      mf.U <- mf.U[rownames( mf.U) %in% ids,, drop=FALSE]
    } else{
      mf.W <- NULL
      mf.U <- stats::model.frame(form3, data = data, na.action = na_rule)
      ids <- c( rownames( mf.U), rownames( mf.X))[duplicated( c( rownames( mf.U), rownames( mf.X)))]  #those rows of data that are good for both parts of the model.
      mf.X <- mf.X[rownames( mf.X) %in% ids,, drop=FALSE]
      mf.U <- mf.U[rownames( mf.U) %in% ids,, drop=FALSE]
    }


  }
  res <- list(ids=ids, mf.X=mf.X, mf.W=mf.W, mf.U=mf.U)

  return( res)
}

"covariate_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))
  # stopifnot(all(is.finite(x)))
}

"get_family_sam" <- function( disty_cases, dist1) {
  error.msg <- paste( c( "Family not implemented. Options are: ", disty_cases, "-- Exitting Now"), collapse=" ")
  disty <- switch( dist1,
                   "bernoulli" = 1,
                   "poisson" = 2,
                   "ippm" = 3,
                   "negative.binomial" = 4,
                   "tweedie" = 5,
                   "gaussian" = 6,
                   "binomial" = 7,
                   {stop( error.msg)} )
  return( disty)
}

"get_logls_sam" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
                            offset, y_is_na, disty, size, powers, control,
                            fits, get_fitted=TRUE){

   if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(y),S))

   #setup the right link function
   if(disty%in%c(1,7)) link <- stats::make.link(link = "logit")
   if(disty%in%c(2,3,4,5)) link <- stats::make.link(link = "log")
   if(disty%in%c(6)) link <- stats::make.link(link = "identity")

   logl_sp <- matrix(NA, nrow=S, ncol=G)

   if(!is.null(U))eta_all <- as.matrix(U) %*% fits$delta
   else eta_all <- rep(0,nrow(X))

    for(ss in 1:S){
      sp_idx<-!y_is_na[,ss]
      for(gg in 1:G){
        eta_mix <- as.matrix(X[sp_idx,]) %*% fits$beta[gg,]
        if(ncol(W)>1) eta_spp <- as.matrix(W[sp_idx,,drop=FALSE]) %*% c(fits$alpha[ss],fits$gamma[ss,])
        else eta_spp <- fits$alpha[ss]
        eta <- eta_spp + eta_mix + eta_all[sp_idx] + offset[sp_idx]
        if(get_fitted) fitted_values[gg,sp_idx,ss] <- link$linkinv(eta)

        if(disty==1) logl_sp[ss,gg] <- sum(dbinom(y[,ss], size =  1, prob =  link$linkinv(eta),log = TRUE))
        if(disty==2) logl_sp[ss,gg] <- sum(dpois(y[,ss], lambda = link$linkinv(eta),log = TRUE))
        if(disty==3) logl_sp[ss,gg] <- y[sp_idx,ss] %*% eta - site_spp_weights[sp_idx,ss] %*% exp(eta)
        if(disty==4) logl_sp[ss,gg] <- sum(dnbinom(y[,ss], mu=link$linkinv(eta), size = exp(-fits$theta[ss]), log=TRUE))
        if(disty==5) logl_sp[ss,gg] <- sum(fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
        if(disty==6) logl_sp[ss,gg] <- sum(dnorm(y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE))
        if(disty==7) logl_sp[ss,gg] <- sum(dbinom(y[,ss],size =  size, prob = link$linkinv(eta),log = TRUE))

      }
      if(!disty%in%3)logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }

   out.list <- list(logl_sp=logl_sp)
   if(get_fitted) out.list$fitted = fitted_values
   return(out.list)
}

"get_logls_spp_sam" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
                            offset, y_is_na, disty, size, powers, control,
                            fits){


  #setup the right link function
  if(disty%in%c(1,7)) link <- stats::make.link(link = "logit")
  if(disty%in%c(2,3,4,5)) link <- stats::make.link(link = "log")
  if(disty%in%c(6)) link <- stats::make.link(link = "identity")

  logl_sp <- matrix(NA, nrow=S, ncol=G)

  if(!is.null(U))eta_all <- as.matrix(U) %*% fits$delta
  else eta_all <- rep(0,nrow(X))

  for(ss in 1:S){
    sp_idx<-!y_is_na[,ss]
    for(gg in 1:G){
      eta_mix <- as.matrix(X[sp_idx,]) %*% fits$beta[gg,]
      if(ncol(W)>1) eta_spp <- as.matrix(W[sp_idx,,drop=FALSE]) %*% c(fits$alpha[ss],fits$gamma[ss,])
      else eta_spp <- fits$alpha[ss]
      eta <- eta_spp + eta_mix + eta_all[sp_idx] + offset[sp_idx]

      if(disty==1) logl_sp[ss,gg] <- sum(dbinom(y[,ss], 1, link$linkinv(eta),log = TRUE))
      if(disty==2) logl_sp[ss,gg] <- sum(dpois(y[,ss], lambda = link$linkinv(eta),log = TRUE))
      if(disty==3) logl_sp[ss,gg] <- y[sp_idx,ss] %*% eta - site_spp_weights[sp_idx,ss] %*% exp(eta)
      if(disty==4) logl_sp[ss,gg] <- sum(dnbinom(y[,ss], mu=link$linkinv(eta), size = exp(-fits$theta[ss]), log=TRUE))
      if(disty==5) logl_sp[ss,gg] <- sum(fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
      if(disty==6) logl_sp[ss,gg] <- sum(dnorm(y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE))
      if(disty==7) logl_sp[ss,gg] <- sum(dbinom(y[,ss], size, link$linkinv(eta),log = TRUE))
    }
    if(!disty%in%3)logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
  }

  ak <- logl_sp + matrix(rep(log(fits$pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  return(sppLogls)

}

"get_offset_sam"  <- function(mf){
  offset <- stats::model.offset(mf)
  if(any(offset!=0))
    return(offset)
  offset <- rep(0, nrow(mf))
  return(offset)
}

"get_power_sam" <- function(disty, power, S){
  if( disty == 5){
    if(length( power) == 1)
      power <- rep(power, S)
    if( length( power) != S)
      stop( "Power parameter(s) not properly specified, exitting now")
  } else {
    power <- -999999
  }
  return( power)
}


"get_site_spp_weights_sam"  <- function(mf, site_spp_weights, sp_names, family){

  # site_spp_wts <- model.weights(mf)
  if(is.null(site_spp_weights))site_spp_weights <- model.weights(mf)

  if(family=='ippm'){
    if(!is.null(site_spp_weights)){
      site_spp_weights <- subset(site_spp_weights, select = colnames(site_spp_weights)%in%sp_names)
    } else {
      site_spp_weights <- matrix(1,nrow(mf),length(sp_names))
    }
  } else {
    if(!is.null(site_spp_weights)){
      site_spp_weights <- replicate(length(sp_names),site_spp_weights)
    } else {
      site_spp_weights <- matrix(1,nrow(mf),length(sp_names))
    }

  }

  return(site_spp_weights)
}

"get_titbits_sam" <- function(titbits, y, X, W, U, spp_weights, site_spp_weights, offset,
                              y_is_na, size, powers, archetype_formula, species_formula, all_formula,
                              control, family)  {
  if( titbits==TRUE)
    titbits <- list(Y = y, X = X, W = W, U = U, spp_weights = spp_weights,
                    site_spp_weights = site_spp_weights, offset = offset,
                    y_is_na = y_is_na, size = size, powers=powers,
                    archetype_formula =  archetype_formula,
                    species_formula = species_formula,
                    all_formula = all_formula,
                    control = control,
                    family = family)
  else{
    titbits <- list()
    if( "Y" %in% titbits)
      titbits$Y <- y
    if( "X" %in% titbits)
      titbits$X <- X
    if( "W" %in% titbits)
      titbits$W <- W
    if( "U" %in% titbits)
      titbits$U <- W
    if( "spp_weights" %in% titbits)
      titbits$spp_weights <- spp_weights
    if( "site_spp_weights" %in% titbits)
      titbits$site_spp_weights <- site_spp_weights
    if( "offset" %in% titbits)
      titbits$offset <- offset
    if( "y_is_na" %in% titbits)
      titbits$y_is_na <- y_is_na
    if( "size" %in% titbits)
      titbits$size <- size
    if( "powers" %in% titbits)
      titbits$powers <- powers
    if( "archetype_formula" %in% titbits)
      titbits$archetype_formula <- archetype_formula
    if( "species_formula" %in% titbits)
      titbits$species_formula <- species_formula
    if( "all_formula" %in% titbits)
      titbits$all_formula <- all_formula
    if( "control" %in% titbits)
      titbits$control <- control
    if( "family" %in% titbits)
      titbits$family <- family
  }
  return( titbits)
}


"get_taus" <- function(pi, logls, G, S){
  fullLogPis <- matrix(rep(log(pi), each=S), nrow=S, ncol=G)
  a_k <- fullLogPis + logls
  a_m <- apply( a_k, 1, max)
  tmp <- exp( a_k - rep( a_m, times=G))
  log_denom <- a_m + log( rowSums( tmp))
  return( exp( a_k - log_denom))
}

"get_X_sam" <- function(archetype_formula, mf.X){
  form.X <- archetype_formula
  form.X[[2]] <- NULL
  form.X <- stats::as.formula(form.X)
  X <- stats::model.matrix(form.X, mf.X)
  tmp.fun <- function(x){ all( x==1)}
  intercepts <- apply( X, 2, tmp.fun)
  X <- X[,!intercepts,drop=FALSE]

  return(as.data.frame(X))
}

"get_W_sam" <- function(species_formula, mf.W){
  form.W <- species_formula
  if(length( form.W)>2)
      form.W[[2]] <- NULL #get rid of outcomes
  W <- as.data.frame(model.matrix( form.W, mf.W))
  colnames(W)[1] <- "const"
  return(W)
}

"get_U_sam" <- function(all_formula, mf.U){
  form.U <- all_formula
  if(!is.null(mf.U)){
    form.W[[2]] <- NULL #get rid of outcomes
  form.U <- update(form.U,~.-1)
  U <- as.data.frame(model.matrix( form.U, mf.U))
  } else {
    U <- NULL
  }
  return(U)
}



"print_input_sam" <- function(y, X, W, U, S, archetype_formula, species_formula,
                              all_formula, family, quiet=FALSE){
  if( quiet)
    return( NULL)
  n.tot <- nrow(y)
  if(family=='ippm'){
    n_pres <- sum(unlist(y)==1,na.rm=TRUE)
    n_bkgrd <- sum(unlist(y[,1])==0,na.rm=TRUE)
    message("There are ", n_pres, " presence observations for ", S," species")
    message("There are ", n_bkgrd, " background (integration) points for each of the ", S," species")
  } else {
    message("There are ", nrow(X), " site observations for ", S," species")
    # message("There are ", ncol(W), " parameters for each species, and ",ncol(X),"parameters for each archetype")
  }

  archetype_formula[[2]] <- NULL
  message("The model for the archetype (grouping) is ", Reduce( "paste", deparse(archetype_formula)))
  if(!is.null(species_formula))
    message("The model for the species is ", Reduce( "paste", deparse(species_formula)))
  if(!is.null(U))
    message("The model for the entire dataset is ", Reduce( "paste", deparse(all_formula)))
  message("You are implementing a ", family, " Species Archetype Model.")
  # else message("You are implementing a ", family, " Partial Species Archetype Model.")
}

"print_starting_values" <-  function(S, G, npx, npw, npu, n, alpha, beta, gamma,
                                     eta, delta, theta){
  message(S, " species.")
  message(G, " groups.")
  message(npx, " archetype coefs.")
  message(npw, " species coefs.")
  message(npu, " all coefs.")
  message(n," sites.")
  message("starting species intercepts:\n",paste(round(alpha,3)," "))
  message("starting archetype parameters:\n",paste(round(beta,3)," "))
  message("starting archetype membership:\n",paste(round(additive_logistic(eta),3)," "))
  message("starting species parameters:\n",paste(round(gamma,3)," "))
  message("starting all parameters:\n", paste(round(delta,3)," "))
  message("starting species specific dispersion parameters:\n",paste(round(theta,3)," "))
}

# "reltol_fun" <- function(logl_n1, logl_n){
#   return(abs(logl_n1 - logl_n) > (abs(logl_n1 - logl_n) / abs(logl_n)))
# }


"sam_random_inits" <- function(alpha, beta, gamma, delta, theta, S, G, X, W, U, disty, mult=0.3, control.sd = control$init_sd){
  if(is.na(control.sd)){
    my.sd <- mult*sd( alpha); if( is.na( my.sd)) my.sd <- 0.1
  } else {
    my.sd <- control.sd
  }
  alpha <- alpha + rnorm(S, sd = my.sd)
  if(is.na(control.sd)){
    my.sd <- mult*sd( beta); if( is.na( my.sd) | my.sd==0) my.sd <- control.sd
  } else {
    my.sd <- control.sd
  }
  beta <- beta + as.numeric(matrix(rnorm(G * ncol(X), mean = 0, sd = my.sd), ncol = ncol(X), nrow = G))
  if( ncol(W)>1){
    if(is.na(control.sd)){
      my.sd <- mult*sd(gamma); if(is.na(my.sd)|my.sd==0) my.sd <- control.sd
    } else {
      my.sd <- control.sd
    }
    gamma <- gamma + as.numeric( matrix(rnorm(S*ncol(W[,-1,drop=FALSE]), mean=0, my.sd), ncol=ncol(W[,-1,drop=FALSE]), nrow=S))
  }
  if( !is.null(U)){
    if(is.na(control.sd)){
      my.sd <- mult*sd(delta); if(is.na(my.sd)|my.sd==0) my.sd <- 0.1
    } else {
      my.sd <- control.sd
    }
    delta <- delta + rnorm(ncol(U), mean=0, my.sd)
  }
  # if(disty %in% c(4,5,6)){
    # my.sd <- mult*sd( theta); if( is.na( my.sd) | my.sd==0) my.sd <- 0.1
    # theta <- theta + as.numeric( rnorm( S, mean=0, my.sd))
  # }
  if(disty %in% c(4,5,6)) return(list(alpha,beta,gamma,delta,theta))
  else return(list(alpha,beta,gamma,delta))
}


"sam_internal_pred_groups" <- function(alpha, beta, taus, gamma, delta,
                                       G, S, X, W, U, offset = NULL, family){

  if (family %in% c("bernoulli","binomial"))
    link.fun <- make.link("logit")
  if (family %in% c("negative.binomial","poisson","tweedie","ippm"))
    link.fun <- make.link("log")
  if (family %in% "gaussian")
    link.fun <- make.link("identity")
  if (is.null(offset))
    offset <- rep(0, nrow(X))

  outpred_arch <- matrix(NA, dim(X)[1], G)
  colnames(outpred_arch) <- paste("G", 1:G, sep = ".")

  if(!is.null(U)) etaAll <- U%*%delta
  else etaAll <- rep(0,nrow(X))

  for (g in seq_len(G)) {
    s.outpred <- matrix(NA, dim(X)[1], length(alpha))
    for (s in seq_len(S)) {
      etaMix <- as.numeric(as.matrix(X)%*%beta[g, ])
      if(ncol(W)>1) etaSpp <- as.numeric(W%*%c(alpha[s],gamma[s, ]))
      else etaSpp <- alpha[s]
      eta <- etaMix + etaSpp + etaAll + offset
      s.outpred[, s] <- link.fun$linkinv(eta)
    }

    outpred_arch[, g] <- apply(s.outpred*rep(taus[, g],each = dim(X)[1]),
                               1, sum)/sum(taus[, g])
  }
  return(outpred_arch)
}

"sam_internal_pred_species" <- function(alpha, beta, taus, gamma, delta,
                                        G, S, X, W, U, offset = NULL, family){

  if (family %in% c("bernoulli","binomial"))
    link.fun <- make.link("logit")
  if (family %in% c("negative.binomial","poisson","ippm"))
    link.fun <- make.link("log")
  if (family %in% "gaussian")
    link.fun <- make.link("identity")
  if (is.null(offset))
    offset <- rep(0, nrow(X))

  outpred_spp <- matrix(0, dim(X)[1], S)

  if(!is.null(U)) etaAll <- U%*%delta
  else etaAll <- rep(0,nrow(X))

  for (g in seq_len(G)) {
    etaMix <- matrix(as.numeric(as.matrix(X)%*%beta[g, ]), nrow(X), S, byrow=FALSE)
    if(ncol(W)>1) etaSpp <- W%*%t(cbind(alpha,gamma))
    else etaSpp <- matrix(alpha, nrow(X), S, byrow=TRUE)
    eta <- etaMix + etaSpp + etaAll + offset
    mug <- link.fun$linkinv(eta)
    outpred_spp <- outpred_spp + mug*matrix(taus[,g], nrow(X), S, byrow=TRUE)
    }

  return(outpred_spp)

}

"setup_inits_sam" <- function(inits, S, G, X, W, U, disty, return_list=TRUE){
  if(is.null(inits))res<-NULL

  npx <- ncol(X)
  if(ncol(W)>1) npw <- ncol(W[,-1])
  else npw <- 0
  if(!is.null(U)) npu <- ncol(U)
  else npu <- 0

  if(is.list(inits)){
    alpha <- as.numeric(inits$alpha)
    beta <- as.numeric(inits$beta)
    eta <- as.numeric(inits$eta)
    if(npw>0){
      gamma <- as.numeric(inits$gamma)
    } else {
      gamma <- rep(-999999,S)
    }
    if(!is.null(U)){
      delta <- inits$delta
    } else{
      delta <- -999999
    }
    if(disty%in%c(4,5,6)){
      theta <- as.numeric(inits$theta)
    } else {
      theta <- rep(-999999,S)
    }
    if(return_list) res <- list(alpha=alpha,beta=beta,eta=eta,gamma=gamma,delta=delta,theta=theta)
    else res <- c(alpha,beta,eta,gamma,delta,theta)
    }

  if(is.numeric(inits)){
    start <- 0
    alpha <- inits[start + 1:S]
    start <- start + S
    beta <- inits[start + 1:((G*npx))]
    start <- start + (G*npx)
    eta <- inits[start + 1:(G - 1)]
    start <- start + (G-1)
    if(npw>0){
      gamma <- inits[start + 1:((S*npw))]
      start <- start + (G*npw)
    } else {
      gamma <- rep(-999999,S)
      # start <- start + S
    }
    if(!is.null(U)){
      delta <- inits[start + 1:npu]
      start <- start + (npu)
    } else{
      delta <- -999999
    }
    if(disty%in%c(4,5,6)){
      theta <- inits[start + 1:S]
    } else {
      theta <- -999999
    }

    if(return_list) res <- list(alpha=alpha,beta=beta,eta=eta,
                                gamma=gamma,delta=delta,theta=theta)
    else res <- c(alpha,beta,eta,gamma,delta,theta)
  }

  return(res)

}

"set_control_sam" <- function(control){
  if (!("maxit" %in% names(control)))
    control$maxit <- 500
  if( !("quiet" %in% names( control)))
    control$quiet <- FALSE
  if( !("cores" %in% names( control)))
    control$cores <- 1
  if (!("trace" %in% names(control)))
    control$trace <- 1
  if( control$quiet)
    control$trace <- 0  #for no tracing
  if (!("init_method" %in% names(control)))
    control$init_method <- 'random2'
  if (!("init_sd" %in% names(control)))
    control$init_sd <- NA
  if (!("minimum_sites_occurrence" %in% names(control)))
    control$minimum_sites_occurrence <- 0
  if (!("ecm_prefit" %in% names(control)))
    control$ecm_prefit <- TRUE
  if (!("ecm_steps" %in% names(control)))
    control$ecm_steps <- 5
  if (!("ecm_refit" %in% names(control)))
    control$ecm_refit <- 1
  if (!("ecm_reltol" %in% names(control)))
    control$ecm_reltol <- 1e-6
  if (!("print_cpp_start_vals" %in% names(control)))
    control$print_cpp_start_vals <- FALSE
  if (!("nreport" %in% names(control)))
    control$nreport <- 10
  if (!("abstol" %in% names(control)))
    control$abstol <- 1e-05
  if (!("reltol" %in% names(control)))
    control$reltol <- sqrt(.Machine$double.eps)
  if (!("optimise_cpp" %in% names( control)))
    control$optimise_cpp <- TRUE
  if (!("getscores_cpp" %in% names( control)))
  control$getscores_cpp <- FALSE
  if (!("loglOnly_cpp" %in% names(control)))
    control$loglOnly_cpp <- FALSE
  if (!("derivOnly_cpp" %in% names( control)))
    control$derivOnly_cpp <- FALSE
  if (!("printparams_cpp" %in% names( control)))
    control$printparams_cpp <- FALSE
  if (!("penalty.pi" %in% names(control)))
    control$penalty.pi <- 0.01
  else
    if (control$penalty.pi < 0) {
      message("Supplied penalty for pis is negative, reverting to the default")
      penalty.pi <- 0.01
    }
  if (!("penalty.alpha" %in% names( control)))
    control$penalty.alpha <- 10
  else
    if (control$penalty.alpha <= 0) {
      message("Supplied penalty for alpha is negative, reverting to the default")
      control$penalty.alpha <- 10
    }
  if( !("penalty.beta" %in% names( control)))
    control$penalty.beta <- 10
  else
    if( control$penalty.beta <=0){
      message("Supplied penalty for betas is negative, reverting to the default")
      control$penalty.beta <- 10
    }
  if( !("penalty.gamma" %in% names( control)))
    control$penalty.gamma <- 10
  else
    if( control$penalty.gamma <=0){
      message("Supplied penalty for gammas is negative, reverting to the default")
      control$penalty.gamma <- 10
    }
  if( !("penalty.delta" %in% names( control)))
    control$penalty.delta <- 10
  else
    if( control$penalty.delta <=0){
      message("Supplied penalty for deltas is negative, reverting to the default")
      control$penalty.delta <- 10
    }
  if( !("penalty.theta" %in% names( control)))
    control$penalty.theta <- c( 10, sqrt( 10))  #the mu and sd of a log-normal
  else
    if( control$penalty.theta[2] <= 0 | length( control$penalty.theta) != 2) {
      message("Supplied penalty parameters for the dispersions is illogical, reverting to the default")
      control$penalty.delta <- c( 10, sqrt( 10))
    }

  return( control)

}

"simulate_ippm_grid" <- function(X,W,n=100,cell_area=1){

              message(paste0("Generating ippm data on a regular ",n,"x",n," grid"))
              x <- y <- 1:n/n
              grid2D <- expand.grid( x, y)
              grid2D$cellArea <- rep(cell_area, nrow( grid2D))  #all cells have same size here
              Xippm <- apply(X,2,function(x)sample(x,n*n,replace = TRUE))
              Wippm <- apply(W,2,function(x)sample(x,n*n,replace = TRUE))
              return(list(grid2D = grid2D, X = Xippm, W = Wippm))
}


"simulate_ippm_outcomes" <- function(X, W, S, grid2D, fitted, cell_area = 1){

  LAMBDAS <- apply(fitted,2,function(x) sum(x * cell_area))
  Ns <- sapply(LAMBDAS, function(x) rpois(n = 1, lambda = x))
  preds_df <- data.frame(idx=1:nrow(X), X)
  presences <- list()

  for(i in seq_len(S)){
    presences[[i]] <- sample(x=preds_df$idx,size=Ns[i], replace=TRUE, prob=fitted[,i]/LAMBDAS[i])
  }

  presence_coords <- lapply(presences,function(x)grid2D[x,1:2])
  presences_sort <- lapply(presences,sort)
  sp_name <- paste0(seq_len(S))
  sp_dat_po_ul<- data.frame(sp=rep(sp_name,unlist(lapply(presences,length))),
                            cell_num=unlist(presences_sort))
  po_matrix <- table_to_species_data(sp_dat_po_ul,
                                     site_id = 'cell_num',species_id = 'sp')
  po_matrix <- po_matrix[,order(as.numeric(as.character(colnames(po_matrix))))]
  po_matrix[po_matrix==0]<-NA
  # po_matrix[order]
  po_covariatesX <- X[as.numeric(rownames(po_matrix)),,drop=FALSE]
  absence_data <- matrix(0,nrow(X),S)
  colnames(po_matrix) <- paste0('spp',sp_name)
  colnames(absence_data) <- paste0('spp',sp_name)
  if(ncol(W)>1){
    po_covariatesW <- W[as.numeric(rownames(po_matrix)),-1,drop=FALSE]
    presence_data <- data.frame(po_matrix,const=1,po_covariatesX,po_covariatesW)
    bkdata <- cbind(absence_data,const=1,X,W[,-1])
  } else {
    presence_data <- data.frame(po_matrix,const=1,po_covariatesX)
    bkdata <- cbind(absence_data,const=1,X)
  }
  # message(colnames(presence_data),colnames(bkdata))
  mm <- rbind(presence_data,bkdata)

  ## calculate out the weights for ippm
  species_specific_cell_counts <- lapply(seq_along(rep(sp_name)),
                                         function(x)table(sp_dat_po_ul[sp_dat_po_ul$sp==rep(sp_name)[x],2]))
  df <- data.frame(id=preds_df$idx,area=grid2D$cellArea)
  sp_weights <- lapply(seq_along(sp_name),
                       function(x)(weights=df$area/as.numeric(species_specific_cell_counts[[x]][match(df$id,as.numeric(names(species_specific_cell_counts[[x]])))])))
  sp_weights_mat <- data.frame(cell_id = 1:nrow(X), do.call(cbind,sp_weights))
  m <- sp_weights_mat
  presence_sites <- m[rowSums(is.na(m[,-1]))!=ncol(m[,-1]), ]
  presence_sites <- data.frame(presence_sites)
  background_sites <- data.frame(cell_id=1:ncol(X),
                                 matrix(rep(grid2D$cellArea,S),nrow(grid2D),S))
  # message(colnames(presence_sites),colnames(background_sites))
  wts <- rbind(presence_sites[,-1],background_sites[,-1])

  return(list(mm=mm,weights=wts))
}

"shrink_taus" <- function(taus, G){
  if( G==1)
    return( taus)
  magical.alpha <- (1-0.8*G)/(0.8*(2-G)-1) ## Dunstan et al., 2013, JABES
  taus_star <- (2*magical.alpha*taus-magical.alpha+1)/(2*magical.alpha - magical.alpha*G + G)
  return(taus_star)
}

"species_mix_boot_parametric" <- function( object, nboot){
  if( nboot > 0){
    if( is.null( object$vcov)){
      message( "An estimate of the variance matrix for regression parameters is required. Please run object$vcov <- vcov(), see ?vcov.species_mix for help")
      return( NULL)
    }
    allCoBoot <- my.rmvnorm( n=nboot, mean=as.numeric(unlist( object$coefs)), sigma=object$vcov, method='eigen')
    return( allCoBoot)
  }
  else{
    boot.estis <- matrix( unlist( object$coef), nrow=1)
    return( boot.estis)
  }
}

"species_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))

  # stopifnot(all(is.finite(x)))
}




"standardise.X" <- function (mat){
  X = scale(as.matrix(mat))
  dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
  dat.sds = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
  return(list(X = X, dat.means = dat.means, dat.sds = dat.sds))
}

"standardise.W" <- function (mat){
  W = scale(as.matrix(mat))
  dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
  dat.sds = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
  return(list(W = W, dat.means = dat.means, dat.sds = dat.sds))
}

"standardise.U" <- function (mat){
  U = scale(as.matrix(mat))
  dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
  dat.sds = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
  return(list(U = U, dat.means = dat.means, dat.sds = dat.sds))
}


"titbit_cleanup" <- function(samfit){

  if(!is.null(samfit$mus))
    samfit$mus <- NULL
  if(!is.null(samfit$scores))
    samfit$scores <- NULL
  if(!is.null(samfit$loglikeSG))
    samfit$loglikeSG <- NULL
  if(!is.null(samfit$loglikeS))
    samfit$loglikeS <- NULL
  if(!is.null(samfit$loglikeSG))
    samfit$loglikeSG <- NULL

  return(samfit)
}


"update_coefs" <- function(old, new, kappa=1){
  if(any(is.na( old)))
    tmp <- new
  else
    tmp <- old + kappa*(new-old)
  return( tmp)
}

