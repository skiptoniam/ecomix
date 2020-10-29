##### Main species mix functions to export #####

#' @title This is how you fit a species archetype model (SAM) in ecomix.
#' @rdname species_mix
#' @name species_mix
#' @description Fits a mixture-of-regressions to identify species archetype
#' models (SAMs).
#' @details species_mix is used to fit mixtures of glms to multivariate
#' species data. The function uses BFGS to optimise the mixture likelihood.
#' There is the option to use EM algorithm to get appropriate starting
#' parameters. `species_mix` acts as a wrapper for species_mix.fit
#' that allows for easier data input. The data frames are merged into
#' the appropriate format for the use in species_mix.fit.
#' Minima is found using vmmin (BFGS). Currently 'bernoulli', 'binomial',
#' 'poisson', 'ippm' (inhomogenous Poisson point process), 'negative.binomial', 'tweedie'
#'  and 'gaussian' distributions can be fitted using the species_mix function.
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
#' @param all_formula an object of class "formula", which is mean to represnet
#'  a constant single set of covariates across all species and groups, typically you might
#'  use this an alternative to an offset, where there might be some bias in the
#'  data which is relatively constant across all species.
#' @param data a matrix of dataframe which contains the 'species_data'
#' matrix, a const and the covariates in the strucute of spp1, spp2, spp3,
#' const, temperature, rainfall. dims of matirx should be
#' nsites*(nspecies+const+covariates).
#' @param nArchetypes The number of archetypes (mixing components/groups) to estimate from the data.
#' @param family The family of statistical family to use within
#' the ecomix models. a  choice between "bernoulli", "poisson", "ippm",
#' "negative.binomial" and "gaussian" familys are possible and applicable
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
#' @param control a list of control parameters for optimisation and calculation.
#' See details. From \code{species_mix.control} for details on optimistaion
#' parameters.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficents. These are family specific, but at a
#' minimum you will need pis (additive_logitic transformed), alpha
#' (intercepts) and beta (mixing coefs).
#' @param standardise Booliean. If TRUE, standarise the covariate data.
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
#' model, "control" for the control arguments used in model fitting, "dist" for
#' the conditional distribution of the species data. Care needs to be taken when
#' using titbits=TRUE in species_mix.multifit(qv) calls as titbits is created
#' for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then
#' setting titbits=TRUE may give users problems with memory.
#' @importFrom graphics abline hist legend lines matplot par plot points polygon
#'  rect
#' @importFrom stats as.formula binomial cooks.distance cov cutree dbinom dist dnbinom dnorm dpois make.link coef glm.fit fitted gaussian glm hclust lm logLik model.matrix model.frame model.offset model.response model.weights pbinom pnbinom pnorm poisson ppois predict qnorm qqnorm quantile rbinom residuals rgamma rnbinom rnorm rpois runif sd uniroot update update.formula nlminb optimise
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
#' simulated_data <- species_mix.simulate(archetype_formula=sam_form,
#' species_formula=sp_form, dat, beta=beta,family="bernoulli")
#' fm1 <- species_mix(sam_form, sp_form, simulated_data,
#'  family = 'bernoulli',  nArchetypes=3)
#'  }

"species_mix" <- function(archetype_formula = NULL,
                          species_formula = stats::as.formula(~1),
                          all_formula = NULL,
                          data, nArchetypes = 3,
                          family="bernoulli", offset=NULL,
                          weights=NULL, bb_weights=NULL, size = NULL, power=1.6,
                          control=NULL, inits=NULL, standardise = FALSE, titbits = TRUE){

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
                         size=size, control=control, inits=inits)

  tmp$dist <- disty_cases[disty]

  if(nArchetypes==1){
    tmp$pis <- tmp$pis
  }else{
    tmp$pis <- additive_logistic(tmp$eta)
  }

  # get logls from parameters

  #calc posterior porbs and pis.
  if(nArchetypes>1){
    fits <- tmp$coefs
    first_fit <- list(y = y, x = X, W = W, U = U,
                      spp_weights = spp_weights,
                      site_spp_weights = site_spp_weights,
                      offset = offset, y_is_na = y_is_na, size = size)
    logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S,
                               disty, get_fitted = FALSE)
    tmp$taus <- get_taus(tmp$pis,logls_mus$logl_sp,G,S)
    tmp$pis <- colSums(tmp$taus)/S
  }

  #Information criteria
  tmp <- calc_info_crit_sam(tmp)

  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_sam(titbits, y, X, W, spp_weights,
                                 site_spp_weights, offset, y_is_na, size,
                                 archetype_formula, species_formula, control,
                                 disty_cases[disty], tmp$removed_species)

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
#'@param y is a matrix genertated from model.response containing the species information. The matrix has the dimensions n_sites * n_species.
#'@param X is a design matrix for the archetype_formula dimension n_sites * n_covariates.
#'@param W is a design matrix for species_formula and will be implemented if species_formula has covariates.
#'@param U is a design matrix fro all_formula and will be implemented if not NULL.
#'@param G is the number of species archetypes that are being estimated.
#'@param S is the number of species to be modelled (this will be calculated internally in species_mix())
#'@param spp_weights These are weights on the species logls and are specifically used in the Bayesian Bootstrap.
#'@param site_spp_weights These are site and species specific weights. For most distributions these will be the same across all species. But this form is required to correctly estiamte the IPPMs. See \link[ecomix]{species_mix} for more details.
#'@param offset this is a vector of site specific offsets, this might be something like area sampled at sites.
#'@param y_is_na This is a logical matrix used specifically with 'ippm' modelling - don't worry about this, it'll be worked out for you. Yay!
#'@param disty the error distribution to used in species_mix estimation. Currently, 'bernoulli', 'poisson', 'ippm' (Poisson point process), 'negative.binomial' and 'guassian' are available - internal conversion of distribution to a integer.
#'@param size The size of each of binomial sample at each site. Length should be the number of sites.
#'@param control this is a list of control parameters that alter the specifics of model fitting. See \link[ecomix]{species_mix.control} for details.
#'@param inits This will be a vector of starting values for species_mix (i.e you've fitted a model and want to refit it).
#'@export

"species_mix.fit" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
                              offset, y_is_na, disty, size, powers, control, inits=NULL){

  if(G==1){
    tmp <- fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights,
                       offset, y_is_na, G, S, disty, size, powers,
                       control=species_mix.control(ecm_refit = 1))
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
                                                control = control)

  } else {
    if(!control$quiet)message('Be careful! You are using your own initial starting values to optimise the species_mix model.')
    inits <- setup_inits_sam(inits, S=S, G=G, X=X, W=W, U=U, disty, return_list = TRUE)
    print(inits)
    starting_values <- inits
  }

  tmp <- sam_optimise(y, X, W, U, offset, spp_weights, site_spp_weights,
                      y_is_na, S, G, disty, size, powers, starting_values, control)

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
#' @param control a list of control parameters for optimisation and calculation.
#' See details. From \code{species_mix.control} for details on optimistaion
#' parameters.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficents. These are distribution specific, but at a
#' minimum you will need pis (additive_logitic transformed), alpha
#' (intercepts) and beta (mixing coefs).
#' @param standardise Booliean. If TRUE, standarise the covariate data.
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
#' model, "control" for the control arguments used in model fitting, "dist" for
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
#' simulated_data <- species_mix.simulate(archetype_formula=sam_form,
#'  species_formula=sp_form, dat,beta=beta,dist="bernoulli")
#' fmods <- species_mix.multifit(sam_form, sp_form, simulated_data,
#'                               family = 'bernoulli', nstart = 10,
#'                               nArchetypes=3)
#' }
"species_mix.multifit" <- function(archetype_formula = NULL,
                                   species_formula = stats::as.formula(~1),
                                   data, nArchetypes = 3, nstart = 10,
                                   mc.cores=1, family="bernoulli",
                                   offset=NULL, weights=NULL, bb_weights=NULL,
                                   size= NULL, control=species_mix.control(),
                                   inits=NULL, standardise = TRUE, titbits = TRUE){

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
  dat <- clean_data_sam(mf, archetype_formula, species_formula, family)

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

  # standarise data if requested.
  x.means <- NULL
  x.sds <- NULL
  w.means <- NULL
  w.sds <- NULL
  if (standardise == TRUE) {
    stand.X <- standardise.X(X)
    X <- as.matrix(stand.X$X)
    if(ncol(W)>1){
      stand.W <- standardise.W(W[,-1,drop=FALSE])
      W <- as.matrix(cbind(1,stand.W$W))
      w.means <- stand.W$dat.means
      w.sds <- stand.W$dat.sds
    }
    x.means <- stand.X$dat.means
    x.sds <- stand.X$dat.sds
  }

  #get family
  disty_cases <- c("bernoulli","poisson","ippm","negative.binomial","tweedie","gaussian","binomial")

  disty <- get_family_sam(disty_cases, family)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # check size
  size <- check_size_binomial(size,nrow(dat$mf.X))

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names,
                                               family)
  spp_weights <- check_spp_weights(bb_weights,S)

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
  if(!control$quiet) print_input_sam(y, X, W, S, archetype_formula,
                                     species_formula, family,
                                     quiet=control$quiet)

  tmp_fun <- function(x){
    if( !control$quiet & nstart>1)
      setTxtProgressBar(pb, x)
        tmpQuiet <- control$quiet
        control$quiet <- TRUE
      # fit species mix.
        tmp <- species_mix.fit(y=y, X=X, W=W, G=nArchetypes, S=S,
                               spp_weights=spp_weights,
                               site_spp_weights=site_spp_weights,
                               offset=offset, disty=disty, y_is_na=y_is_na,
                               size=size, control=control, inits=inits)

        tmp$dist <- disty_cases[disty]

        if(nArchetypes==1) tmp$pis <- tmp$pis
        else tmp$pis <- additive_logistic(tmp$eta)

        #calc posterior porbs and pis.
        if(nArchetypes>1)
          tmp$taus <- calc_post_probs_sam(tmp$pis,tmp$loglikeSG)

        tmp$pis <- colSums(tmp$taus)/S

        #Information criteria
        tmp <- calc_info_crit_sam(tmp)

        #titbits object, if wanted/needed.
        tmp$titbits <- get_titbits_sam(titbits, y, X, W, spp_weights,
                                       site_spp_weights, offset, y_is_na, size,
                                       archetype_formula, species_formula, control,
                                       disty_cases[disty], tmp$removed_species)
      class(tmp) <- c("species_mix")
      return( tmp)
   }

  #    require( parallel)
  if( !control$quiet & nstart>1)
    pb <- txtProgressBar(min = 1, max = nstart, style = 3, char = "c[_] ")

   #Fit the model many times
   many_starts <- surveillance::plapply(seq_len(nstart), tmp_fun,
                                        .parallel = mc.cores,
                                        .verbose = !control$quiet)

   class(many_starts) <- c("species_mix.multifit")

   return(many_starts)
}

#'@title species_mix.control
#'@rdname species_mix.control
#'@name species_mix.control
#'@description This is the control function used to tweak the species_mix model.
#'@param quiet Should any reporting be performed? Default is FALSE, for reporting.
#'@param cores The number of cores to use in fitting of species_mix models. These will be largely used to model the species-specific parameteres.
#'@param init_method The method to use for initialisation. The options are "random2", "kmeans", "kmed". The default uses random2, which is a kmeans with noise added to the cluster.
#'@param init_sd The amount of noise to add to the initailisation the default is 1.
#'@param minimum_sites_occurrence a integer which determins the number of minimum sites present for any species for it to be included in the initialisation step. This removes rare species from initial groupings. They are then included in the overall analysis.
#'@param ecm_prefit Logical if TRUE the model will run a slower EM algorithm fit to find starting values.
#'@param ecm_steps int Default is 3, the number of EM iterations to get to starting values.
#'@param ecm_refit int Default is 1, number of times to refit using EM.
#'@param ecm_reltol A function or value which gives the tolerance in the EM loglikeihood estimation.
#'@param theta_range Two positive values use as penalities for estimating the dispersion parameters (theta) in a negative.binomial SAM or PSAM.
#'@param pen_param A penality for the em fitting.
#'@param update_kappa Penalities for how fast parameters update during each EM step.
#'@param print_cpp_start_vals A call to check what parameter estimates are being passed to C++ for optimisation.
#'@param maxit_cpp The number of iterations to run in C++. Default is 1000.
#'@param trace_cpp Non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information.
#'@param nreport_cpp The number of iterations to report the loglikelihood. Default is 10.
#'@param abstol_cpp Absolute tolerance. Defaults to 0 so the absolute convergence test is not used. If the objective function is known to be non-negative, the previous default of 1e-20 would be more appropriate.
#'@param reltol_cpp Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.
#'@param conv_cpp Has the model convered previously.
#'@param printparams_cpp Print the parameter estimates within C++.
#'@param optimise_cpp Should optimisation for estimation occur? If TRUE (default) optimisation will occur. If FALSE no optimisation is performed.
#'@param loglOnly_cpp Should the log-likelihood be calculated? If TRUE (default) then log-likelihood is calculated and returned. If FALSE then the log-likelihood is not calculated for return.
#'@param derivOnly_cpp Should the scores be evaluated at the (final) parameter values. If TRUE (default) then they are calculated. If FALSE then they are not calculated.
#'@param getscores_cpp Return scores.
#'@param \dots Other control calls.

#'@export
"species_mix.control" <- function(quiet = FALSE,
                                  cores = 1,
                                  ## intialisation controls
                                  init_method = 'random2',
                                  init_sd = NULL,
                                  minimum_sites_occurrence = 0,
                                  init_glmnet = FALSE,
                                  ## EM algorithim controls
                                  ecm_prefit = TRUE,
                                  ecm_steps = 5,
                                  ecm_refit = 3,
                                  ecm_reltol = 1e-6,
                                  ## partial mixture penalities
                                  theta_range = c(0.001,10),
                                  pen_param = 1.25,
                                  update_kappa = c(1,0.5,1),
                                  ## c++ controls
                                  print_cpp_start_vals = FALSE,
                                  maxit_cpp = 1000,
                                  trace_cpp = 1,
                                  nreport_cpp = 10,
                                  abstol_cpp = sqrt(.Machine$double.eps),
                                  reltol_cpp = sqrt(.Machine$double.eps),
                                  conv_cpp = 1,
                                  printparams_cpp = 0,
                                  optimise_cpp = 1,
                                  loglOnly_cpp = 0,
                                  derivOnly_cpp = 0,
                                  getscores_cpp = 0,
                                  ...){
               #general controls
  rval <- list(quiet = quiet, cores = cores,
               #initialisation controls
               init_method = init_method, init_sd = init_sd,
               minimum_sites_occurrence = minimum_sites_occurrence,
               init_glmnet = init_glmnet,
               #em controls
               ecm_prefit = ecm_prefit,
               ecm_refit = ecm_refit,
               ecm_steps = ecm_steps,
               ecm_reltol = ecm_reltol,
               # partial controls
               theta_range = theta_range,
               pen_param = pen_param,
               update_kappa = update_kappa,
               #cpp controls
               print_cpp_start_vals = print_cpp_start_vals,
               maxit_cpp = maxit_cpp, trace_cpp = trace_cpp,
               nreport_cpp = nreport_cpp,
               abstol_cpp = abstol_cpp, reltol_cpp = reltol_cpp,
               conv_cpp = conv_cpp,
               printparams_cpp = printparams_cpp, optimise_cpp = optimise_cpp,
               loglOnly_cpp = loglOnly_cpp, derivOnly_cpp = derivOnly_cpp,
               getscores_cpp = getscores_cpp)
  rval <- c(rval, list(...))
  if (is.null(rval$ecm_reltol))
    rval$ecm_reltol <- sqrt(.Machine$double.eps)
  rval
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

  if(object$dist=='ippm'){
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
    disty <- get_family_sam(disty_cases, object$dist)
    dumbOut <- capture.output(
      samp.object <- species_mix.fit(y=object$titbits$Y,
                                     X=object$titbits$X,
                                     W=object$titbits$W,
                                     offset = object$titbits$offset,
                                     spp_weights = all.wts[dummy,,drop=TRUE],
                                     site_spp_weights = object$titbits$site_spp_weights,
                                     G = object$G,
                                     S = object$S,
                                     y_is_na = object$titbits$y_is_na,
                                     disty = disty,
                                     size = object$titbits$size,
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
    tmp <- surveillance::plapply(seq_len(nboot), my.fun, .parallel = mc.cores)
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
#' @param dat a matrix of variables to simulate data from.
#' @param nArchetypes number of groups to simulate.
#' @param alpha coefficents for each species archetype. vector S long.
#' @param beta coefficents for each species archetype. Matrix of G x number of
#'  parameters. Each row is a different species archetype.
#' @param gamma coefficents for each species archetype. Matrix of S x number of
#'  parameters. Each row is a different species archetype.
#' @param theta coefficents for the dispersion variables for negative.binomial
#' and gaussian distributions - should be number of species long
#' @param size Is for the binomial model and this respresents the number of binomial trials per site, can be fixed or vary.
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
#'                                        dat, nArchetypes = 4, beta=beta,
#'                                        family="bernoulli")
#' }

## need to update this to take the new formula framework and simulate ippm data.
"species_mix.simulate" <-  function(archetype_formula,
                                    species_formula,
                                    all_formula=NULL,
                                    dat,
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
    U <- stats::model.matrix(all_formula, dat)
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

  X <- stats::model.matrix(archetype_formula, dat)
  W <- stats::model.matrix(species_formula, dat)
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

##### S3 class SAM functions #####
#' @rdname AIC.species_mix
#' @name AIC.species_mix
#' @title Return AIC from a species_mix model
#' @param object A species mix object
#' @param k AIC parameter
#' @param \\dots Ignored
#' @export

"AIC.species_mix" <- function (object, k=NULL, ...){

  effect.param <- length(object$beta) + G + S + ifelse(disty %in%  "negative.binomial",S,0)

  p <- length(unlist(object$coefs))
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


## This needs fixing because it includes the null coefs.
"BIC.species_mix" <-  function (object, ...){
  p <- length(unlist(object$coefs))
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
  if(object$dist%in%c('negative.binomial','gaussian')){
    res$theta <- object$coef$theta
    names(res$theta) <- object$names$spp
  }
  return(res)
}

#' @rdname plot.species_mix
#' @name plot.species_mix
#' @title plot.species_mix
#' @param x a fitted species_mix model.
#' @param species which species residuals to plot. Default is "AllSpecies".
#' @param fitted.scale log or logit, this enables the plotting of residuals to be on the linear predictor scale.
#' @param \\dots Extra plotting arguments.
#' @details Plot random quantile residuals (RQR). "RQR" produces residuals for each species.
#' @references  Dunn, P.K. and Smyth G.K. (1996) Randomized Quantile Residuals. Journal of Computational and Graphical Statistics \emph{5}: 236--244.
#' @export

"plot.species_mix" <- function (x,
                                species="AllSpecies",
                                fitted.scale="response",
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
    preds <- sam_internal_pred_species(x$coef$alpha, x$coef$beta, x$taus,
                              x$coef$gamma, x$G, x$S, x$titbits$X,
                              x$titbits$W, x$titbits$offset, x$dist)

    preds <- preds[,sppID]

    switch( fitted.scale,
            log = { loggy <- "x"},
            logit = { loggy <- ""; preds <- log( preds / (1-preds))},
            {loggy <- ""})
    plot( preds, obs.resid, xlab="Fitted", ylab="RQR", main="Residual versus Fitted", sub="Colours separate species", pch=20, col=rep( 1:S, each=x$n), log=loggy)
    abline( h=0)

  # }
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
#'@param prediction_type Do you want to produce 'archetype' or 'species' level predictions. default is 'archetype'.
#'@param \\dots Ignored
#'@description Predict species archetypes from a species_mix model. You can also predict the conditional species predictions using "prediction_type='species'".
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
#' simulated_data <- species_mix.simulate(archetype_formula=sam_form,
#'  species_formula=sp_form, dat,beta=beta,dist="bernoulli")
#' fm1 <- species_mix(sam_form, sp_form, simulated_data,
#'  family = 'bernoulli',  nArchetypes=3)
#'preds_fm1 <- predict(fm1)
#'}

"predict.species_mix" <- function (object, object2 = NULL, newdata = NULL,
                                   offset = NULL, nboot = 0, alpha = 0.95,
                                   mc.cores = 1, prediction_type='archetype', ...){
  if (is.null(newdata)) {
    X <- object$titbits$X
    W <- object$titbits$W
    offset <- object$titbits$offset
  } else {
    arch.fm <- as.formula(object$titbits$archetype_formula)
    spp.fm <- as.formula(object$titbits$species_formula)
    if (length(arch.fm) == 3) arch.fm[[2]] <- NULL
    arch.fm <- update(arch.fm,~.-1)
    X <- model.matrix(arch.fm, as.data.frame(newdata))
    W <- model.matrix(spp.fm, as.data.frame(newdata))
    offset <- model.frame(arch.fm, data = newdata)
    offset <- model.offset(offset)
  }

  if (is.null(offset))
    offset <- rep(0, nrow(X))

  S <- object$S
  G <- object$G
  n <- nrow(X)
  npx <- object$npx
  npw <- object$npw

  spp_wts <- object$titbits$spp_weights
  site_spp_wts <- object$titbits$site_spp_weights

  disty_cases <- c("bernoulli","poisson","ippm","negative.binomial","tweedie","gaussian","binomial")
  disty <- get_family_sam(disty_cases, object$titbits$family)
  taus <- object$taus
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

  alphaIn <- c(NA, as.numeric(object$coefs$alpha))
  alphaIn <- alphaIn[-1]
  betaIn <- c(NA, as.numeric(object$coef$beta))
  betaIn <- betaIn[-1]
  # etaIn <- c(NA, as.numeric(object$coef$eta))
  # etaIn <- etaIn[-1]
  if (npw>0) {
    gammaIn <- c(NA, as.numeric(object$coef$gamma))
    gammaIn <- gammaIn[-1]
    usetheta <- 1
  } else {
    gammaIn <- -999999
    usetheta <- 0
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
      tmp <- switch (prediction_type,
                     archetype = sam_internal_pred_groups(alpha = object$coefs$alpha,
                               beta = object$coefs$beta,
                               gamma = object$coefs$gamma,
                               taus = taus, G = G, S = S, X = X, W = W,
                               offset = offset, family = object$dist),
                     species = sam_internal_pred_species(alpha = object$coefs$alpha,
                                                        beta = object$coefs$beta,
                                                        gamma = object$coefs$gamma,
                                                        taus = taus, G = G, S = S, X = X, W = W,
                                                        offset = offset, family = object$dist))
    } else {
      nboot <- segments[seg]
      bootSampsToUse <- (sum( segments[1:seg])-segments[seg]+1):sum(segments[1:seg])

      # add in species level preds.
      tmp <- lapply(bootSampsToUse,function(ii)switch (prediction_type,
                                                       archetype = sam_internal_pred_groups(alpha = alphaBoot[ii,],
                                                                 beta = matrix(betaBoot[ii,],G,npx),
                                                                 gamma = matrix(gammaBoot[ii,],S,npw),
                                                                 taus = taus, G = G, S = S, X = X, W = W,
                                                                 offset = offset, family = object$dist),
                                                       species = sam_internal_pred_species(alpha = alphaBoot[ii,],
                                                                                           beta = matrix(betaBoot[ii,],G,npx),
                                                                                           gamma = matrix(gammaBoot[ii,],S,npw),
                                                                                           taus = taus, G = G, S = S, X = X, W = W,
                                                                                           offset = offset, family = object$dist)))

    }

    if (nboot == 0) {
      ret_grp <- tmp
      if(prediction_type%in%"archetype")colnames(ret_grp) <- object$names$SAMs
      if(prediction_type%in%"species")colnames(ret_grp) <- object$names$spp
      return(ret_grp)
    }

    if(prediction_type%in%"archetype"){
      bootPreds <- matrix(do.call("cbind",lapply(tmp,c)), nrow = nrow(X) * G,  ncol = nboot)
    }
    if(prediction_type%in%"species"){
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
    if(prediction_type%in%"archetype") tmp <- matrix(row.exp, nrow = nrow(X), ncol = G)
    if(prediction_type%in%"species") tmp <- matrix(row.exp, nrow = nrow(X), ncol = S)
    bPreds$fit <- tmp
    tmp.grp <- sweep(bootPreds, 1, row.exp, "-")
    tmp.grp <- tmp.grp^2
    tmp.grp <- sqrt(rowSums(tmp.grp)/(nboot - 1))
    if(prediction_type%in%"archetype") tmp.grp <- matrix(tmp.grp, nrow = nrow(X), ncol = G)
    if(prediction_type%in%"species") tmp.grp <- matrix(tmp.grp, nrow = nrow(X), ncol = S)
    bPreds$ses <- tmp.grp
    if(prediction_type%in%"archetype") colnames(bPreds$fit) <- colnames(bPreds$ses) <- object$names$SAMs
    if(prediction_type%in%"species") colnames(bPreds$fit) <- colnames(bPreds$ses) <- object$names$spp
    tmp.fun <- function(x) return(quantile(bootPreds[x, ],
                                           probs = c(0, alpha) + (1 - alpha)/2,
                                           na.rm = TRUE))
    tmp1 <- parallel::mclapply(seq_len(nrow(bootPreds)), tmp.fun,
                               mc.cores = mc.cores)
    tmp1 <- do.call("rbind", tmp1)
    if(prediction_type%in%"archetype"){
      tmp1 <- array(tmp1, c(nrow(X), G, 2), dimnames = list(NULL,
                                                          NULL, NULL))
      bPreds$cis <- tmp1[, 1:G, ]
      dimnames(bPreds$cis) <- list(NULL, object$names$SAMs, c("lower", "upper"))
    }
    if(prediction_type%in%"species"){
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
    if(prediction_type%in%"archetype"){
      tmp1 <- array(tmp1, c(nrow(X), G, 2), dimnames = list(NULL,
                                                          NULL, NULL))
    bPreds$fit_cis <- tmp1[, 1:G, ]
    dimnames(bPreds$fit_cis) <- list(NULL, object$names$SAMs, c("lower", "upper"))
    }
    if(prediction_type%in%"species"){
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
#' simulated_data <- species_mix.simulate(archetype_formula=sam_form,
#'  species_formula=sp_form, dat,beta=beta,dist="bernoulli")
#' fm1 <- species_mix(sam_form, sp_form, simulated_data,
#'  family = 'bernoulli',  nArchetypes=3)
#'print(fm1)}

"print.species_mix" <-  function (x,...){
  cat(x$titbits$family, "species_mix model\n")
  cat("\nMixing probabilities\n")
  print(x$pi)
  cat("\nCoefficents\n")
  print(x$coef)

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
  cat(x[[1]]$titbits$family,"species_mix model\n")
  cat("\nMixing probabilities\n")
  print(x$pi)
  cat("\nCoefficents\n")
  print(x$coef)

}

#' @title Estimate residuals for a species_mix object
#' @rdname residuals.species_mix
#' @name residuals.species_mix
#' @param object A returned species_mix model object.
#' @param \dots additional calls for residual function
#' @param type The type of residuals to estimate. Default is "RQR" (Random Quantile Residuals).
#' But you can also simulate many Random Quantile Residuals using "SimRQR".
#' @param control Default uses species_mix.control()
#' @export
#' @description  The randomised quantile residuals ("RQR", from Dunn and Smyth, 1996) are defined by
#' their marginal distribution function (marginality is over other species observations within that site;
#' see Woolley et al, in prep).

"residuals.species_mix" <- function( object, ..., type="RQR", control=species_mix.control()) {
    if( ! type %in% c("RQR","deviance","pearson"))
      stop( "Unknown type of residual requested.\n Only deviance and RQR (for randomised quantile residuals) are implemented\n")
    if(type=="pearson"){
      stop("pearson residuals not implemented yet.")
    }
    if( type=="RQR"){
      resids <- matrix( NA, nrow=object$n, ncol=object$S)
      switch( object$dist,
              bernoulli = { fn <- function(y,mu,logtheta) pbinom( q=y, size=1, prob=mu, lower.tail=TRUE)},
              poisson = { fn <- function(y,mu,logtheta) ppois( q=y, lambda=mu, lower.tail=TRUE)},
              ippm = { fn <- function(y,mu,logtheta) ppois( q=y, lambda=mu, lower.tail=TRUE)},
              negative.binomial = { fn <- function(y,mu,logtheta) pnbinom( q=y, mu=mu, size=1/exp( logtheta), lower.tail=TRUE)},
              gaussian = { fn <- function(y,mu,logtheta) pnorm( q=y, mean=mu, sd=exp( logtheta), lower.tail=TRUE)})


      for( ss in 1:object$S){
        if( object$dist %in% c("bernoulli","poisson","ippm","negative.binomial")){
          tmpLower <- fn( object$titbits$Y[,ss]-1, object$mus[,ss,], object$coef$theta[ss])
          tmpUpper <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$theta[ss])
          tmpLower <- rowSums( tmpLower * object$pis)
          tmpLower <- ifelse( tmpLower<0, 0, tmpLower) #get rid of numerical errors for really small negative values
          tmpLower <- ifelse( tmpLower>1, 1, tmpLower) #get rid of numerical errors for 1+epsilon.
          tmpUpper <- rowSums( tmpUpper * object$pis)
          tmpUpper <- ifelse( tmpUpper<0, 0, tmpUpper) #get rid of numerical errors for really small negative values
          tmpUpper <- ifelse( tmpUpper>1, 1, tmpUpper) #get rid of numerical errors for 1+epsilon.
          resids[,ss] <- runif( object$n, min=tmpLower, max=tmpUpper)
          resids[,ss] <- qnorm( resids[,ss])
        }
        if( object$dist == "gaussian"){
          tmp <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$theta[ss])
          tmp <- rowSums( tmp * object$pis)
          resids[,ss] <- qnorm( tmp)
        }
      }
      if( !control$quiet & sum( resids==Inf | resids==-Inf | is.na(resids))>0)
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


#'@rdname vcov.species_mix
#'@name vcov.species_mix
#'@title Estimate the Variance-covariance matrix for a species_mix object.
#'@description Calculates variance-covariance matrix from a species_mix object
#'@param object an object obtained from fitting a RCP (for region of common profile) mixture model. Such as that generated from a call to species_mix(qv).
#'@param object2 an object of class \code{species_mix} containing bootstrap samples of the parameter estimates (see species_mix_boot(qv)). If NULL (default) the bootstrapping is performed from within the vcov function. If not null, then the vcov estimate is obtained from these bootstrap samples.
#'@param method the method to calculate the variance-covariance matrix. Options are:'FiniteDifference' (default), \code{BayesBoot}, \code{SimpleBoot}, and \code{EmpiricalInfo}. The two bootstrap methods (\code{BayesBoot} and \code{SimpleBoot}, see species_mix_boot(qv)) should be more general and may possibly be more robust. The \code{EmpiricalInfo} method implements an empirical estimate of the Fisher information matrix, I can not recommend it however. It seems to behave poorly, even in well behaved simulations. It is computationally thrifty though.
#'@param nboot the number of bootstrap samples to take for the bootstrap estimation. Argument is ignored if !method \%in\% c(\code{FiniteDifference},'EmpiricalInfo').
#'@param mc.cores the number of cores to distrbute the calculations on. Default is 4. Set to 1 if the computer is running Windows (as it cannot handle forking -- see mclapply(qv)). Ignored if method=='EmpiricalInfo'.
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
#' simulated_data <- species_mix.simulate(archetype_formula=sam_form,
#'  species_formula=sp_form, dat,beta=beta,dist="bernoulli")
#' fm1 <- species_mix(sam_form, sp_form, simulated_data,
#'  family = 'bernoulli',  nArchetypes=3)
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
    # if( Sys.info()['sysname'] == "Windows")
    X <- object$titbits$X
    W <- object$titbits$W
    offset <- object$titbits$offset
    spp_wts <- object$titbits$spp_weights
    site_spp_wts <- object$titbits$site_spp_weights
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
    control <- object$titbits$control

    # values for optimisation.
    inits <- object$coefs
    start_vals <- setup_inits_sam(inits, S, G, X, W, disty, return_list = TRUE)

    # parameters to optimise
    alpha <- as.numeric(start_vals$alpha)
    beta <- as.numeric(start_vals$beta)
    eta <- as.numeric(start_vals$eta)
    gamma <- as.numeric(start_vals$gamma)
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
      gamma.score <- -999999
    }
    if(disty%in%c(4,6)){
      control$optiDisp <- as.integer(1)
      theta.score <- as.numeric(rep(NA, length(theta)))
    }else{
      control$optiDisp <- as.integer(0)
      theta.score <- -999999
    }
    scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,gamma.score,theta.score))))
    conv <- FALSE

    #model quantities
    pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
    mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
    loglikeS <- as.numeric(rep(NA, S))
    loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

    #remove finite for now, as we only need it for ippm
    if (method %in% c("FiniteDifference")) {
      grad_fun <- function(x) {
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
          gamma <- -999999
          # start <- start + 1
        }
        if(disty%in%c(4,6)){
          theta <- x[start + seq_len(S)]
        } else {
          theta <- -999999
        }
        #c++ call to optimise the model (needs pretty good starting values)
        tmp <- .Call("species_mix_cpp",
                     as.numeric(as.matrix(y)), as.numeric(as.matrix(X)), as.numeric(as.matrix(W[,-1,drop=FALSE])), as.numeric(offset), as.numeric(spp_wts),
                     as.numeric(as.matrix(site_spp_wts)), as.integer(as.matrix(!y_is_na)), as.numeric(size),
                     # SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_weights, SEXP Rsite_spp_weights, SEXP Ry_not_na, // data
                     as.integer(S), as.integer(G), as.integer(npx), as.integer(npw), as.integer(n),
                     as.integer(disty),as.integer(control$optiDisp),as.integer(control$optiPart),
                     # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
                     as.double(alpha), as.double(beta), as.double(eta), as.double(gamma), as.double(theta),
                     # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
                     alpha.score, beta.score, eta.score, gamma.score, theta.score, as.integer(getscores), scores,
                     # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
                     pis_out, mus, loglikeS, loglikeSG,
                     # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
                     as.integer(control$maxit_cpp), as.integer(control$trace_cpp), as.integer(control$nreport_cpp),
                     as.numeric(control$abstol_cpp), as.numeric(control$reltol_cpp),  as.integer(conv), as.integer(control$printparams_cpp),
                     # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
                     as.integer(0), as.integer(0), as.integer(1),
                     # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
                     PACKAGE = "ecomix")

        tmp1 <- c(alpha.score, beta.score, eta.score)
        if( npw > 0)#class( object$titbits$species_formula) == "formula")
          tmp1 <- c(tmp1, gamma.score)
        if(disty%in%c(4,6))
          tmp1 <- c( tmp1, theta.score)
        return(tmp1)
      }
      hess <- numDeriv::jacobian(grad_fun, x=unlist(object$coefs),method="simple")
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

###### ECM for fitting SAMs ######
"fit.ecm.sam" <- function(y, X, W, U=NULL, spp_weights, site_spp_weights, offset,
                          y_is_na, G, S, disty, size, powers, control, starting.sam = NULL){

  n <- nrow(y)
  starting.sam <-  NULL
  bestOfAllMods <- list( logl=-Inf)
  for(t in seq_len(control$ecm_refit)) {
    message(paste0("ECM restart ",t," of ", control$ecm_refit,""))
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
                                           control = species_mix.control(init_method = "random2"))

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
      fits$alpha <- ecomix:::update_coefs(fits$alpha,new.sp.params[,1])
      if(ncol(W)>1){
        fits$gamma <- ecomix:::update_coefs(fits$gamma,new.sp.params[,-1])
      } else {
        fits$gamma <- rep(-99999,S)
      }

      ## update the dispersion parameters.
      if(disty %in% c(4,5,6)) {
        fm_theta <- plapply(seq_len(S), apply_optimise_thetaParams,
                            y, X, W, U, taus, fits,
                            site_spp_weights, offset,
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
        new.deltas <- optimize(f = llogl.allParams, interval = c(-30,30), maximum = TRUE,
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
                             fits, powers),
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
      mod <- list(logl = new.logl,alpha = fits$alpha, beta = fits$beta, gamma = fits$gamma,
                  delta = fits$delta, theta = fits$theta, pis = fits$pis, taus = taus)
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
    }
    if( mod$logl > bestOfAllMods$logl)
      bestOfAllMods <- mod
  }

  return(bestOfAllMods)
}

"e.step" <- function(y, X, W, U, site_spp_weights, offset,
                     y_is_na, disty,#data
                     fits, sp.powers,
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
        out.taus[ss,gg] <- (sum(dbinom(y[,ss], 1, p = check.p, log = TRUE)))
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
        out.taus[ss,gg] <- (sum(dbinom(y[,ss], size, p = check.p, log = TRUE)))
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
  update.sppParams <- optimize(f = llogl.sppParams, interval = c(-30,30), maximum = TRUE, ss = ss,
                               y = y, X = X, W = W, U=U, taus = taus, fits = fits,
                               site_spp_weights = site_spp_weights,
                               offset = offset, y_is_na = y_is_na,
                               disty = disty, size = size, powers = powers)$maximum
  return(update.sppParams)
}

"apply_optimise_thetaParams" <- function(ss, y, X, W, U, taus, fits,
                                   site_spp_weights, offy, y_is_na,
                                   disty, size, powers){
  update.theta <- optimize(f = llogl.thetaParams, interval = c(-30,30), maximum = TRUE, ss = ss,
                           y = y, X = X, W = W, U=U, taus = taus, fits = fits,
                           site_spp_weights = site_spp_weights,
                           offy = offset, disty = disty, size = size, powers = powers)$maximum
  return(update.theta)
}

"apply_optimise_betas" <- function(gg, y, X, W, U, site_spp_weights, offset,
                                   y_is_na, disty, taus, fits, size, powers){

  n <- nrow(y)
  if(disty %in% c(2,3,4,5))
    glmnet.family <- "poisson"; glm.family <- poisson();
    if(disty %in% c(1,7))
      glmnet.family <- "binomial"; glm.family <- binomial();
      if(disty %in% c(6))
        glmnet.family <- "gaussian"; glm.family <- gaussian();

        if(disty %in% 4) get.mus <- e.step(y, X, W, U, site_spp_weights,
                                           offset, y_is_na, disty,
                                           fits, powers,  get.fitted = TRUE)$fitted

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
          Y_s <- as.vector(Y_s/site_weights)
        }
        if(disty%in%c(7)){
          Y_s <- as.matrix(cbind(Y_s,size_s-Y_s))
        }

        if(ncol(X_s)==1) X_s <- cbind(1,X_s)

        fit1 <- glmnet::glmnet(x = as.matrix(X_s), y = Y_s, family = glmnet.family, weights = obs.weights+1e-6, offset = offy, nlambda = 100, intercept = FALSE)
        new.betas <- coef(fit1)[,ncol(coef(fit1))][-1]
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
  if(disty%in%6) link <- make.link('identity')

  if(!is.null(U)) all.etas <- as.matrix(U)%*%c(fits$delta)
  else all.etas <- rep(0,n)

  mix.etas <- as.matrix(X)%*%t(fits$betas)

  for(gg in seq_len(G)) {
    cw.eta <- as.matrix(W)%*%x + mix.etas[,gg] + all.etas + offset
    if(disty %in%  1)
      out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = 1, p =  link$linkinv(cw.eta), log = TRUE))
    if(disty %in%  2)
      out <- out + sum(taus[ss,gg]*dpois(y[,ss], lambda =  link$linkinv(cw.eta), log = TRUE))
    if(disty %in%  3)
      out <- out + sum(taus[ss,gg]*(y[sp_idx,ss] %*% cw.eta - site_spp_weights[sp_idx,ss] %*%  link$linkinv(cw.eta)))
    if(disty %in%  4)
      out <- out + sum(taus[ss,gg]*dnbinom(y[,ss], mu =  link$linkinv(cw.eta), size = 1/fits$theta[ss], log = TRUE))
    if(disty %in%  5)
      out <- out + sum(taus[ss,gg]*fishMod::dTweedie(y[,ss], mu =  link$linkinv(cw.eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
    if(disty %in%  6)
      out <- out + sum(taus[ss,gg]*dnorm(y[,ss], mean =  link$linkinv(cw.eta), sd = sqrt(fits$theta[ss]), log = TRUE))
    if(disty %in%  7)
      out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = size, p = link$linkinv(cw.eta), log = TRUE))
  }
  return(out)
}


"llogl.thetaParams" <- function(x, ss, y, X, W, U, taus, fits, site_spp_weights, offset, disty, size, powers) {

  n <- nrow(y)
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
     out <- out + sum(taus[ss,gg]*fishMod::dTweedie(y[,ss], mu = link$linkinv(eta), phi = x, p = powers[ss], LOG = TRUE));
   if(disty==6)
     out <- out + sum(taus[ss,gg]*dnorm(y[,ss], mean = link$linkinv(eta), sd = sqrt(x), log = TRUE))
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
      spp.etas <- W %*% c(fits$alpha[ss],fits$gamma[ss,])
    }else{
      spp.etas <- W %*% c(fits$alpha[ss])
    }

    for(gg in seq_len(G)) {
      cw.eta <- spp.etas + mix.etas[,gg] + all.etas + offset
      if(disty %in%  1)
        out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = 1, p =  link$linkinv(cw.eta), log = TRUE))
      if(disty %in%  2)
        out <- out + sum(taus[ss,gg]*dpois(y[,ss], lambda =  link$linkinv(cw.eta), log = TRUE))
      if(disty %in%  3)
        out <- out + sum(taus[ss,gg]*(y[sp_idx,ss] %*% cw.eta - site_spp_weights[sp_idx,ss] %*%  link$linkinv(cw.eta)))
      if(disty %in%  4)
        out <- out + sum(taus[ss,gg]*dnbinom(y[,ss], mu =  link$linkinv(cw.eta), size = 1/fits$theta[ss], log = TRUE))
      if(disty %in%  5)
        out <- out + sum(taus[ss,gg]*fishMod::dTweedie(y[,ss], mu =  link$linkinv(cw.eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
      if(disty %in%  6)
        out <- out + sum(taus[ss,gg]*dnorm(y[,ss], mean =  link$linkinv(cw.eta), sd = sqrt(fits$theta[ss]), log = TRUE))
      if(disty %in%  7)
        out <- out + sum(taus[ss,gg]*dbinom(y[,ss], size = size, p = link$linkinv(cw.eta), log = TRUE))
    }

  }
  return(out)
}



###### SAM internal functions for fitting ######
"starting_values_wrapper" <- function(y, X, W, U, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, size, control){
  if(isTRUE(control$ecm_prefit)){
    if(!control$quiet)message('Using ECM algorithm to find starting values; using ',
                              control$ecm_refit,'refits')
    emfits <- fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights,
                            offset, y_is_na, G, S, disty, size, powers, control)
    # bf <- which.max(vapply(emfits,function(x)c(x$logl),c(logl=0)))
    # emfit <- emfits[[bf]]
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
                                              # spp_weights = spp_weights,
                                              site_spp_weights = site_spp_weights,
                                              offset = offset, y_is_na = y_is_na,
                                              G = G, S = S,
                                              disty=disty,size=size,powers = powers,
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
  return(start_vals)
}




# "apply_optimise_spp_theta" <- function(ss, first_fit, fits,
#                                        G, disty, pis,
#                                        theta.range = c(0.0001,100)){
#   thet <- optimise(f = theta.logl, interval = theta.range, ss, first_fit,
#                    fits, G, disty, pis, theta.range,
#                    maximum = TRUE)$maximum
#   return(thet)
# }
#
# "theta.logl" <- function(theta, ss, first_fit, fits, G,
#                          disty, pis, theta.range){
#
#   if(disty==4)link <- make.link('log')
#   if(disty==6)link <- make.link('identity')
#
#   if(ncol(first_fit$W)>1){
#     eta_spp <- first_fit$W %*% c(fits$alpha[ss],fits$gamma[ss,])
#   }else{
#     eta_spp <- first_fit$W %*% c(fits$alpha[ss])
#   }
#
#   if(!is.null(first_fit$U)){
#     eta_all <- first_fit$U %*% fits$delta
#   } else {
#     eta_all <- rep(0,nrow(first_fit$y))
#   }
#
#   eta_mix <- first_fit$x %*% t(fits$beta)
#   offy <- first_fit$offset
#   y <- first_fit$y[,ss]
#   logls <- rep( NA, G)
#   for( gg in 1:G){
#     eta <- eta_spp + eta_mix[,gg] + eta_all + offy
#     mus <- link$linkinv(eta)
#     if(disty==4)logls[gg] <- sum( dnbinom(x = y, mu = mus, size = theta, log=TRUE))
#     if(disty==6)logls[gg] <- sum( dnorm(x = y, mean = mus, sd = theta, log=TRUE))
#   }
#   ak <- logls + log(pis)
#   am <- max(ak)
#   ak <- exp( ak-am)
#   sppLogls <- am + log( sum( ak))
#
#   return(sppLogls)
# }

# "incomplete_negbin_logl" <- function( x, pis, first_fit, fits, G, S) {
#   fits$beta <- matrix(x, nrow=nrow( fits$beta), ncol=ncol( fits$beta))
#   tmp <- get_incomplete_logl_nb( pis, first_fit, fits, G, S)
#   return( -tmp)
# }
#
#
# "get_incomplete_logl_nb" <- function(pis, first_fit, fits, G, S, theta.range=c(0.001, 10), pen.parm=1.25){
#   logls <- matrix(NA, nrow=S, ncol=G)
#   if(ncol(first_fit$W)>1){
#     eta_spp <- first_fit$W %*% t(cbind(fits$alpha,fits$gamma))
#   } else {
#     eta_spp <- first_fit$W %*% t(fits$alpha)
#   }
#   eta_mix <- first_fit$x %*% t(fits$beta)
#   incomplete.logl <- 0
#   for( ss in seq_len(S)){
#     for( gg in seq_len(G)){
#       eta <- eta_spp[,ss] + eta_mix[,gg] + first_fit$offset
#       logls[ss,gg] <- sum(dnbinom(first_fit$y[,ss], mu=exp(eta), size=exp(-fits$theta[ss]), log=TRUE))
#     }
#   }
#   ak <- logls + matrix( rep( log( pis), each=S), nrow=S, ncol=G)
#   am <- apply( ak, 1, max)
#   ak <- exp( ak-am)
#   sppLogls <- am + log( rowSums( ak))
#   # sppLogls <- sppLogls  + dbeta((exp(-fits$theta)-theta.range[1]) / (theta.range[2]- theta.range[1]), pen.parm, pen.parm, log=TRUE)
#   logl <- sum( sppLogls)
#   return( logl)
# }

## function for starting values using penalities
"get_initial_values_sam" <- function(y, X, W, U=NULL, site_spp_weights,
                                   offset, y_is_na, G, S, disty, size, powers, control) {

  n <- nrow(y)
  prev.min <- floor(n*control$minimum_sites_occurrence);
  sel.omit.spp <- which(colSums(y>0) <= prev.min)
  if(length(sel.omit.spp)==0) sel.omit.spp <- -1*(1:S)

  starting.sam <- list(alpha = rep(0,S), theta = rep(1,S));
  message("Initialising starting values")


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
    fit1 <- many.fitt(y, X, W, U, site_spp_weights,
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
    fit1 <- mvabund::manylm(tmpform,data = datX)
    starting.sam$theta <- (deviance(fit1)/fit1$df.residual)
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
    mrwdist <- kmed::distNumeric(spp.beta, sp.beta, method = "mrw")
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
    random_coefs <- ecomix:::sam_random_inits(alpha = starting.sam$alpha,beta = grp_coefs,
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
  starting.sam$taus <- ecomix:::shrink_taus(taus, G=G)
  starting.sam$pis <- colMeans(taus)

  return(starting.sam)
}

"many.fit" <- function(y, X, W, U, site_spp_weights, offset, y_is_na, G, S, disty, size, powers, control){

  options(warn = -1)
  fm_sp_mods <-  ecomix:::plapply(seq_len(S), apply_glmnet_sam_inits, y, X, W, U,
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


"apply_glmnet_sam_inits" <- function(ss, y, X, W, U = NULL, site_spp_weights,
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
  if (disty==7){
    outcomes <- as.matrix(cbind(y[ids_i,ss],size[ids_i]-y[ids_i,ss]))
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

  if(!disty %in% c(5)){

    lambda.seq <- sort(unique( c( seq( from=1/0.001, to=1, length=25),
                                   seq( from=1, to=.1, length=10))),
                       decreasing=TRUE)
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
  } else { #Tweedie needs an unconstrained fit.  May cause problems in some cases, especially if there is quasi-separation...
      df3 <- data.frame(y=outcomes, offy=offset, Intercept= 1, df)
      tmp.fm1 <- fishMod::tglm(y~-1+.-offy+offset( offy),
                               wts=c(site_spp_weights[,ss]),
                               data=df3, p=power[ss], vcov=FALSE,
                               residuals=FALSE, trace=0)
      my_coefs <- t(as.matrix(tmp.fm1$coef,ncol=1))
      theta <- log(tmp.fm1$coef["phi"])
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

# "apply_glm_sam_inits" <- function(ss, y, X, W, U=NULL, site_spp_weights,
#                                   offset, y_is_na, disty, size){
#
#   # which family to use?
#   if(disty == 1 | disty == 7 ) #binomials
#     fam <- binomial() #glmnet
#   if(disty == 2 | disty == 3 | disty == 4)
#     fam <- poisson()#"poisson"
#   if(disty == 6)
#     fam <- gaussian()#"gaussian"
#
#   #set up the index of NAs (this is for ppms)
#   ids_i <- !y_is_na[,ss]
#
#   if (disty==3){
#     outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
#   } else {
#     outcomes <- as.matrix(y[ids_i,ss])
#   }
#   if (disty==7){
#     outcomes <- as.matrix(cbind(y[ids_i,ss],size[ids_i]-y[ids_i,ss]))
#   }
#
#   mm <- cbind(W[ids_i,-1,drop=FALSE],X[ids_i,,drop=FALSE])
#   if(!is.null(U)) mm <- cbind(mm,U[ids_i,,drop=FALSE])
#
#   offy <- as.numeric(offset[ids_i])
#   wts <- as.numeric(site_spp_weights[ids_i,ss])
#   dat <- data.frame(outcomes=outcomes, mm, offy=offy)
#
#   if( disty %in% c(1,2,3,4,6,7)){
#     tmpform <- as.formula(paste0('outcomes ~ 1 +', paste0(colnames(mm),collapse = "+"),' + offset(offy)'))
#     ft_sp <- try(glm(formula = tmpform, data = dat, weights=wts, family=fam), silent=FALSE)
#     if (any(class(ft_sp)[1] %in% 'try-error')){
#       my_coefs <- rep(NA, ncol(X[ids_i,]))
#       names(my_coefs) <- colnames(W[ids_i,,drop=FALSE],X[ids_i,,drop=FALSE],U)
#     } else {
#       # if(ncol(X)==1) my_coefs <- t(as.matrix(my.coefs[-1]))
#       my_coefs <- coef(ft_sp)
#     }
#
#     ##estimate the starting dispersion parameter.
#     theta <- -99999
#     if( disty == 4){
#       tmp <- MASS::theta.mm(outcomes,ft_sp$fitted.values,
#                             weights=as.matrix(site_spp_weights[ids_i,ss]),
#                             dfr=length(outcomes), eps=1e-4)
#       if(tmp>2) tmp <- 2
#       theta <- tmp#log( 1/tmp)
#     }
#     if( disty == 6){
#       preds <- as.numeric(ft_sp$linear.predictors)
#       theta <- log( sqrt( sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard
#     }
#   }
#   # species intercpets
#   alpha <- my_coefs[1]
#   # mixture coefs
#   beta <- my_coefs[match(colnames(X), names(my_coefs))]
#   # species coefs apart from intercept
#   if(ncol(W)>1) gamma <-  my_coefs[match(colnames(W[,-1,drop=FALSE]), names(my_coefs))]
#   else gamma <- -99999
#
#   if(!is.null(U)) delta <-  my_coefs[match(colnames(U), colnames(my_coefs))]
#   else delta <- -99999
#
#
#   return(list(alpha = alpha, beta = beta, gamma = gamma, delta = delta, theta = theta))
# }

# #get the conditional maxima for species specific parameters.
# "apply_glm_spp_coefs_sams" <- function(ss, y, X, W, U=NULL,
#                                        G, taus,
#                                        site_spp_weights,
#                                        offset, y_is_na, disty, fits, size){
#   if(disty %in% c(1,7))
#     fam <- binomial()
#   if(disty %in% c(2,3,4))
#     fam <- poisson()
#   if(disty %in% c(6))
#     fam <- gaussian()
#
#   ids_i <- !y_is_na[,ss]
#
#   if (disty %in% c(7)){
#     outcomes <- as.matrix(cbind(y[ids_i,ss],size[ids_i]-y[ids_i,ss]))
#   }
#   if (disty %in% c(2,4)){
#     outcomes <-y[ids_i,ss]
#   }
#   if (disty == 3){
#     outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
#
#   }
#   if (disty == 6){
#     outcomes <- y[ids_i,ss, drop=FALSE]
#   }
#
#   out1 <- kronecker(rep( 1, G),  outcomes)
#   W1 <- kronecker(rep( 1, G), W[ids_i,,drop=FALSE])
#   wts1 <- kronecker(rep( 1, G), as.numeric(site_spp_weights[ids_i,ss]))*rep(taus[ss,],each=length(site_spp_weights[ids_i,ss]))
#   offyArea <- kronecker(rep( 1, G), offset[ids_i])
#   offyGrp <- X[ids_i,] %*% t(fits$beta)
#   offyGrp <- as.numeric(offyGrp)
#
#   if(!is.null(U)){
#     U1 <- kronecker(rep( 1, G), U[ids_i,,drop=FALSE])
#     offyAll <- U1 %*% (fits$delta)
#   } else {
#     offyAll <- rep(0,nrow(W1))
#   }
#   offyAll <- as.numeric(offyAll)
#   offy <- offy1 + offyGrp + offyAll
#
#   # length(offy1)
#   # length(offyAll)
#   # dim(out1)
#   # dim(W1)
#   # dim(wts1)
#   if(disty %in% c(1,2,3,6,7)){
#     tmpform <- as.formula('out1 ~ -1+W1+offset( offyArea)+offset( offyGrp)+offset( offyAll)')
#     ft_sp <- try(stats::glm( tmpform, weights=wts1, family=fam),silent = TRUE)
#     # ft_sp <- try(stats::glm.fit(x=as.data.frame(W1),
#     #                             y=as.vector(out1),
#     #                             weights=as.numeric(wts1),
#     #                             offset=as.numeric(offy), #start = fits$alpha[ss],
#     #                             family=fam), silent=FALSE)
#     if(!ft_sp$converged){
#       tmpform <- as.formula('out1 ~ -1+W1+offset( offy1)+offset( offyGrp)+offset( offyAll)')
#       ft_sp <- try(mgcv::gam( tmpform, weights=wts1, family=fam), silent=FALSE)
#
#     }
#     if (class(ft_sp)[1] %in% 'try-error'){
#        my_coefs <- rep(NA, ncol(W1))
#     } else {
#       my_coefs <- stats::coef(ft_sp)
#       names(my_coefs) <- c(colnames(y)[ss],colnames(W[,-1,drop=FALSE]))
#     }
#   }
#   if(disty %in% c(4)){
#     tmpform <- as.formula( paste('out1','-1+W1+offset(offy)', sep='~'))
#     ft_sp <- try(mgcv::gam(tmpform, weights=wts1, family=mgcv::negbin(theta=exp(-fits$theta[ss]))))
#     kount1 <- 1
#     while( any(class( ft_sp)[1] %in% 'try-error') & kount1 < 10){
#       kount1 <- kount1 + 1
#       theta <- 10 * exp(-fits$theta[ss])
#       ft_sp <- try(mgcv::gam(tmpform, weights=wts1, family=mgcv::negbin(theta=theta)))
#       }
#     my_coefs <- ft_sp$coefficients
#     names(my_coefs) <- c(colnames(y)[ss],colnames(W[,-1,drop=FALSE]))
#   }
#   return(list(alpha = my_coefs[1], gamma = my_coefs[-1]))
# }
#
# #get the conditional maximums for group coefs.
# "apply_glm_mix_coefs_sams" <- function(gg, y, X, W, U = NULL,
#                                        site_spp_weights, offset,
#                                        y_is_na, disty,
#                                        taus, fits, mus, size){
#
#   ### setup the data stucture for this model.
#   # dim(y)
#   Y_s <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
#   size_s <- matrix(rep(size,ncol(y)),nrow(y),ncol(y))[!y_is_na]
#   X_no_NA <- list()
#   for (jj in 1:ncol(y)){
#     X_no_NA[[jj]] <- X[!y_is_na[,jj],,drop=FALSE]
#   }
#   X_s <- do.call(rbind, X_no_NA)
#
#   if(!is.null(U)){
#     U_no_NA <- list()
#     for (jj in 1:ncol(y)){
#       U_no_NA[[jj]] <- U[!y_is_na[,jj],,drop=FALSE]
#     }
#   U_s <- do.call(rbind, U_no_NA)
#   offy3 <- U_s %*% fits$delta
#   } else {
#   offy3 <- rep(0,nrow(X_s))
#   }
#
#   n_ys <- sapply(X_no_NA,nrow)
#
#
#   ## set up the weights
#   tau.weights <- rep(taus[,gg,drop=FALSE],c(n_ys))
#   if(disty %in% 4) tau.weights <- rep(taus[,gg,drop=FALSE],c(n_ys))/(1+rep(fits$theta,n_ys)*as.vector(mus[k,,]))
#   site.weights <- as.matrix(as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
#   obs.weights <- tau.weights*site.weights
#   # obs.weights <- obs.weights+1e-6
#
#   ## set up the offsets.
#   offy_mat <- replicate(ncol(y),offset)
#
#   offy1 <- as.matrix(unlist(as.data.frame(offy_mat[!y_is_na])))
#   if(ncol(W)>1){
#     offy2 <- W %*% t(cbind(fits$alpha,fits$gamma))
#   } else {
#     offy2 <- W %*% c(fits$alpha)
#   }
#   offy2 <- as.matrix(unlist(as.data.frame(offy2[!y_is_na])))
#   offy <- as.numeric(offy1 + offy2 + offy3)
#
#   # which family to use?
#   if( disty == 1 | disty == 7)
#     fam <- binomial()
#   if( disty == 2 | disty == 3 |disty == 4)
#     fam <- poisson()
#   if( disty == 6)
#     fam <- gaussian()
#   if (disty==3){
#     Y_s <- as.matrix(Y_s/site_weights)
#   } else {
#     Y_s <- as.matrix(Y_s)
#   }
#   if(disty%in%c(7)){
#     Y_s <- as.matrix(cbind(Y_s,size_s-Y_s))
#   }
#
#   if(disty %in% c(1,2,3,6,7)){ #don't use for tweedie - try and fit negative_binomial using glm.fit.nbinom
#     # ft_mix <- try(glm(Y_s ~ X_s - 1, offset = offy, weights = wts_sXsite_weights, family = fam))
#     ft_mix <- try(glm.fit(x = as.data.frame(X_s),
#                           y =Y_s,
#                           weights = c(wts_sXsite_weights),
#                           offset = offy,
#                           family = fam), silent = TRUE)
#     if (class(ft_mix)[1] %in% 'try-error'){
#       mix_coefs <- rep(NA, ncol(X_s))
#     } else {
#       mix_coefs <- ft_mix$coefficients#[-1]
#     }
#   }
#
#   return(c(mix_coefs))
# }

# # just need one glmo
# "glm_all_coefs_sams" <- function(y, X, W, U, site_spp_weights, offset,
#                                  y_is_na, disty, taus, fits, mus, size){
#
#   ### setup the data stucture for this model.
#   ### This is going to be a huge model over S*G
#   G <- ncol(taus)
#   S <- nrow(taus)
#   Y_allspp <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
#   size_allspp <- as.matrix(matrix(rep(size,ncol(y)),nrow(y),ncol(y))[!y_is_na], ncol=1)
#   X_no_NA <- list()
#   U_no_NA <- list()
#   for (jj in 1:ncol(y)){
#     X_no_NA[[jj]] <- X[!y_is_na[,jj],,drop=FALSE]
#     U_no_NA[[jj]] <- U[!y_is_na[,jj],,drop=FALSE]
#   }
#   X_allspp <- do.call(rbind, X_no_NA)
#   U_allspp <- do.call(rbind, U_no_NA)
#   n_ys <- sapply(X_no_NA,nrow)
#
#   ## Set up the giant glm
#   Y_SG <- kronecker(rep(1,G),Y_allspp)
#   size_SG <- kronecker(rep(1,G),size_allspp)
#   U_SG <- kronecker(rep(1,G),U_allspp)
#   tausSG <- rep(as.vector(taus),rep(c(n_ys),G))
#
#   site_weights <- kronecker(rep(1,G),as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
#   wts_tausXsite_weights <- tausSG*site_weights  ## weird things could be happening here ippms
#   offy_mat <- replicate(ncol(y)*G,offset)
#   offy1 <- unlist(as.data.frame(offy_mat[!y_is_na]))
#
#
#   ## species etas
#   if(ncol(W)>1) offy2 <- W %*% t(cbind(fits$alpha,fits$gamma))
#   else offy2 <- W %*% c(fits$alpha)
#   offy2 <- unlist(as.data.frame(offy2[!y_is_na]))
#   offy2 <- kronecker(rep(1,G),offy2)
#
#   ## mix etas
#   offy3 <- X_allspp %*% t(fits$beta)
#   offy3 <- unlist(as.data.frame(offy2[!y_is_na]))
#
#   # sum all the etas together.
#   offy <- as.numeric(offy1 + offy2 + offy3)
#
#
#   # which family to use?
#   if( disty == 1 | disty == 7)
#     fam <- binomial()
#   if( disty == 2 | disty == 3 |disty == 4)
#     fam <- poisson()
#   if( disty == 6)
#     fam <- gaussian()
#   if (disty==3){
#     Y_SG <- as.matrix(Y_SG/site_weights)
#   } else {
#     Y_SG <- as.matrix(Y_SG)
#   }
#   if(disty%in%c(7)){
#     Y_SG <- as.matrix(cbind(Y_SG,size_SG-Y_SG))
#   }
#
#   if(disty %in% c(1,2,3,6,7)){ #don't use for tweedie - try and fit negative_binomial using glm.fit.nbinom
#     ft_all <- try(glm.fit(x = as.data.frame(U_SG),
#                           y = Y_SG,
#                           weights = c(wts_tausXsite_weights),
#                           offset = offy,
#                           family = fam), silent = TRUE)
#     if (class(ft_all)[1] %in% 'try-error'){
#       all_coefs <- rep(NA, ncol(U))
#     } else {
#       all_coefs <- ft_all$coefficients
#     }
#   }
#   return(c(mix_coefs))
# }


# "get_incomplete_logl_sam_V2" <- function(eta, first_fit, fits,
#                                       spp_weights, G, S, disty){
#
#   pis <- additive_logistic(eta)
#   # if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(first_fit$y),S))
#
#   #setup the right link function
#   if(disty%in%c(1,7)) link <- stats::make.link(link = "logit")
#   if(disty%in%c(2,3,4,5)) link <- stats::make.link(link = "log")
#   if(disty%in%c(6)) link <- stats::make.link(link = "identity")
#
#   logl_sp <- matrix(NA, nrow=S, ncol=G)
#
#   if(!is.null(first_fit$U))eta_all <- as.matrix(first_fit$U) %*% fits$delta
#   else eta_all <- rep(nrow(first_fit$x),0)
#
#   for(ss in 1:S){
#     sp_idx<-!first_fit$y_is_na[,ss]
#     for(gg in 1:G){
#       eta_mix <- as.matrix(first_fit$x[sp_idx,]) %*% fits$beta[gg,]
#       if(ncol(first_fit$W)>1) eta_spp <- as.matrix(first_fit$W[sp_idx,,drop=FALSE]) %*% c(fits$alpha[ss],fits$gamma[ss,])
#       else eta_spp <- fits$alpha[ss]
#       eta <- eta_spp + eta_mix + eta_all[sp_idx] + first_fit$offset[sp_idx]
#       # if(get_fitted) fitted_values[gg,sp_idx,ss] <- link$linkinv(eta)
#       if(disty==1) logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta),log = TRUE))
#       if(disty==2) logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss], lambda = link$linkinv(eta),log = TRUE))
#       if(disty==3) logl_sp[ss,gg] <- first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% exp(eta)
#       if(disty==4) logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss], mu=link$linkinv(eta), size = exp(-fits$theta[ss]), log=TRUE))
#       if(disty==5) stop('No tweedie')
#       if(disty==6) logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE))
#       if(disty==7) logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], first_fit$size, link$linkinv(eta),log = TRUE))
#     }
#     if(!disty%in%3)logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
#   }
#
#   ak <- logl_sp + matrix(rep(log( pis), each=S), nrow=S, ncol=G)
#   am <- apply( ak, 1, max)
#   ak <- exp( ak-am)
#   sppLogls <- am + log( rowSums( ak))
#   # theta.range=c(0.001, 10); pen.parm=1.25
#   # if(disty==4)  sppLogls <- sppLogls  + dbeta((exp(-fits$theta[ss])-theta.range[1]) / (theta.range[2]- theta.range[1]), pen.parm, pen.parm, log=TRUE)
#   logl <- sum( sppLogls)
#   return(logl)
# }

# "get_initial_values_sam" <- function(y, X, W, U, spp_weights, site_spp_weights,
#                                      offset, y_is_na, G, S, disty, size, control){
#
#   # get intial model fits
#   starting_values <- initiate_fit_sam(y, X, W, U, spp_weights, site_spp_weights,
#                                       offset, y_is_na, G, S, disty, size, control)
#
#   #if any are errors then remove them from the models forever.
#   # updated_y <- update_species_data_structure(y, y_is_na,
#   #                                            spp_weights,
#   #                                            site_spp_weights,
#   #                                            starting_values$species_to_remove)
#   # y <- updated_y[[1]]
#   # y_is_na <- updated_y[[2]]
#   # spp_weights <- updated_y[[3]]
#   # site_spp_weights <- updated_y[[4]]
#   # S <- ncol(y)
#
#   fits <- list(alpha=starting_values$alpha,
#                beta=starting_values$beta,
#                gamma=starting_values$gamma,
#                delta=starting_values$delta,
#                theta=starting_values$theta)
#   first_fit <- list(y = y, x = X, W = W, U = U, spp_weights = spp_weights,
#                     site_spp_weights = site_spp_weights, offset = offset,
#                     y_is_na = y_is_na, size = size)
#
#   res <- list()
#   res$fits <- fits
#   res$first_fit <- first_fit
#   res$taus <- starting_values$taus
#   res$pis <- colSums(starting_values$taus)/S
#   return(res)
# }

# "get_logls_sam" <- function(first_fit, fits, spp_weights, G, S, disty, get_fitted=TRUE){
#
#    if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(first_fit$y),S))
#
#    #setup the right link function
#    if(disty%in%c(1,7)) link <- stats::make.link(link = "logit")
#    if(disty%in%c(2,3,4,5)) link <- stats::make.link(link = "log")
#    if(disty%in%c(6)) link <- stats::make.link(link = "identity")
#
#    logl_sp <- matrix(NA, nrow=S, ncol=G)
#
#    if(!is.null(first_fit$U))eta_all <- as.matrix(first_fit$U) %*% fits$delta
#    else eta_all <- rep(nrow(first_fit$x),0)
#
#     for(ss in 1:S){
#       sp_idx<-!first_fit$y_is_na[,ss]
#       for(gg in 1:G){
#         eta_mix <- as.matrix(first_fit$x[sp_idx,]) %*% fits$beta[gg,]
#         if(ncol(first_fit$W)>1) eta_spp <- as.matrix(first_fit$W[sp_idx,,drop=FALSE]) %*% c(fits$alpha[ss],fits$gamma[ss,])
#         else eta_spp <- fits$alpha[ss]
#         eta <- eta_spp + eta_mix + eta_all[sp_idx] + first_fit$offset[sp_idx]
#         if(get_fitted) fitted_values[gg,sp_idx,ss] <- link$linkinv(eta)
#         # logl_sp[ss,gg] <- switch(disty,
#         #                          1 = sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta),log = TRUE)),
#         #                          2 = sum(dpois(first_fit$y[,ss], lambda = link$linkinv(eta),log = TRUE)),
#         #                          3 = first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% exp(eta),
#         #                          4 =  sum(dnbinom(first_fit$y[,ss], mu=link$linkinv(eta), size = exp(-fits$theta[ss]), log=TRUE)),
#         #                          5 = stop('No tweedie'),
#         #                          6 = sum(dnorm(first_fit$y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE)),
#         #                          7 = sum(dbinom(first_fit$y[,ss], first_fit$size, link$linkinv(eta),log = TRUE)))
#
#         if(disty==1) logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(eta),log = TRUE))
#         if(disty==2) logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss], lambda = link$linkinv(eta),log = TRUE))
#         if(disty==3) logl_sp[ss,gg] <- first_fit$y[sp_idx,ss] %*% eta - first_fit$site_spp_weights[sp_idx,ss] %*% exp(eta)
#         if(disty==4) logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss], mu=link$linkinv(eta), size = exp(-fits$theta[ss]), log=TRUE))
#         if(disty==5) stop('No tweedie')
#         if(disty==6) logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=eta,sd=exp(fits$theta[ss]),log=TRUE))
#         if(disty==7) logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], first_fit$size, link$linkinv(eta),log = TRUE))
#       }
#       if(!disty%in%3)logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
#     }
#
#    out.list <- list(logl_sp=logl_sp)
#    if(get_fitted) out.list$fitted = fitted_values
#    return(out.list)
# }

# "fitmix_ECM_sam" <- function(y, X, W, U, spp_weights, site_spp_weights, offset,
#                              y_is_na, G, S, disty, size, control){
#
#   ite <- 1
#   restart_ite <- 1
#   logl_old <- -99999999
#   logl_new <- -88888888
#
#   # get starting values
#   starting_values <- get_initial_values_sam(y = y, X = X, W = W, U = U,
#                                             spp_weights = spp_weights,
#                                             site_spp_weights = site_spp_weights,
#                                             offset = offset, y_is_na = y_is_na,
#                                             G = G, S = S,
#                                             disty = disty, size = size,
#                                             control = control)
#
#   # first e-step
#   fits <- starting_values$fits
#   taus <- starting_values$taus
#   if(!disty==3) taus <- ecomix:::shrink_taus(taus, max_tau = 1/G + 0.1, G)
#   pis <- starting_values$pis
#   first_fit <- starting_values$first_fit
#   logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
#
#   while(control$ecm_reltol(logl_new,logl_old) & ite <= control$ecm_steps){
#     if(restart_ite>10){
#       message('cannot find good starting values with initialisation\n
#               and random starting values\n\n
#               Please check the number of groups and coefs.')
#       break
#     }
#
#     pis <- colSums(taus)/S
#
#     if (any(pis < sqrt(.Machine$double.eps))) {
#        if(restart_ite==1){
#          message('Pis have gone to zero - restarting with new initialisation \n')
#          starting_values <- get_initial_values_sam(y = y, X = X, W = W, U=U,
#                                                 spp_weights = spp_weights,
#                                                 site_spp_weights = site_spp_weights,
#                                                 offset = offset, y_is_na = y_is_na,
#                                                 G = G, S = S,
#                                                 disty = disty, size=size,
#                                                 control = control)
#       pis <- starting_values$pis
#       fits <- starting_values$fits
#       taus <- starting_values$taus
#       first_fit <- starting_values$first_fit
#       } else {
#       message('Pis have gone to zero - restarting with random inits \n')
#       taus <- matrix(runif(S*G),S); taus <- taus/rowSums(taus);
#       pis <- colSums(taus)/S;
#       fits$beta <- matrix(rnorm(G*(ncol(X))),G,(ncol(X)))
#       if(ncol(W)>1)fits$gamma <- matrix(rnorm(S*(ncol(W[,-1,drop=FALSE]))),S,ncol(W[,-1,drop=FALSE]))
#       if(!is.null(U))fits$delta <- rnorm(ncol(U))
#       fits$alpha <- rnorm(S)
#       if (disty%in%c(4,6)) fits$theta <- rep(0.05,S)
#       }
#       ite <- 1
#       restart_ite <- restart_ite + 1
#     }
#
#     # m-step
#     ## update the betas
#     if(disty==4){
#       fm_mix_coefs <- nlminb(start=fits$beta, objective=incomplete_negbin_logl, gradient=NULL,
#                              hessian=NULL, pis=pis, first_fit=first_fit, fits=fits, G=G, S=S)
#       fm_mix_coefs_mat <- matrix(fm_mix_coefs$par,G,ncol(X))
#     } else {
#       fm_mix_coefs <- surveillance::plapply(seq_len(G),
#                                             apply_glm_mix_coefs_sams,
#                                             y, X, W, site_spp_weights,
#                                             offset, y_is_na, disty, taus,
#                                             fits, logls_mus$fitted, size,
#                                             .parallel = control$cores,
#                                             .verbose = FALSE)
#
#       fm_mix_coefs_mat <- do.call(rbind,fm_mix_coefs)
#       # print(fm_mix_coefs_mat)
#     }
#     fits$beta <- update_coefs(fits$beta, fm_mix_coefs_mat)
#
#     ## check this one out.
#     fm_spp_coefs <- surveillance::plapply(seq_len(S),
#                                           apply_glm_spp_coefs_sams,
#                                           y, X, W, U, G, taus, site_spp_weights,
#                                           offset, y_is_na, disty, fits, size,
#                                           .parallel = control$cores,
#                                           .verbose = FALSE)
#
#     # get and update the intercepts.
#     alpha <- unlist(lapply(fm_spp_coefs, `[[`, 1))
#     fits$alpha <- update_coefs(fits$alpha,alpha)
#
#     # get and update the gamma if included in the model.
#     if(ncol(W)>1) {
#       gamma <- do.call(rbind,lapply(fm_spp_coefs, `[[`, 2))
#       fits$gamma <- update_coefs(fits$gamma,gamma)
#     } else {
#       fits$gamma <- -99999
#     }
#
#     ## update the betas
#     if(disty==4){
#       fm_mix_coefs <- nlminb(start=fits$beta, objective=incomplete_negbin_logl, gradient=NULL,
#                     hessian=NULL, pis=pis, first_fit=first_fit, fits=fits, G=G, S=S)
#       fm_mix_coefs_mat <- matrix(fm_mix_coefs$par,G,ncol(X))
#     } else {
#       fm_mix_coefs <- surveillance::plapply(seq_len(G),
#                                           apply_glm_mix_coefs_sams,
#                                           y, X, W, U,
#                                           site_spp_weights,
#                                           offset, y_is_na, disty, taus,
#                                           fits, logls_mus$fitted, size,
#                                           .parallel = control$cores,
#                                           .verbose = FALSE)
#
#       fm_mix_coefs_mat <- do.call(rbind,fm_mix_coefs)
#     }
#     fits$beta <- update_coefs(fits$beta, fm_mix_coefs_mat)
#
#     ## update the deltas
#     if(!is.null(U)){
#       fm_all_coefs <- glm_all_coefs_sams(y, X, W, U,
#                                        site_spp_weights,
#                                        offset, y_is_na, disty, taus,
#                                        fits, logls_mus$fitted, size)
#       fits$delta <- update_coefs(fits$delta, fm_all_coefs)
#
#     }
#
#
#     ## need a function here that updates the dispersion parameter.
#     if(ite >= init_steps){
#       if(disty%in%c(4)){
#         fits$theta <- exp(-fits$theta)
#         fm_theta <- surveillance::plapply(seq_len(S), apply_optimise_spp_theta,
#                                          first_fit, fits,
#                                          G, disty, pis,
#                                          .parallel = control$cores,
#                                          .verbose = FALSE)
#         theta <- unlist(lapply(fm_theta, `[[`, 1))
#         theta <- log(1/theta)
#         fits$theta <- update_coefs(log(1/fits$theta),theta,control$update_kappa[2])
#       }
#
#       if(disty%in%c(6)){
#         fits$theta <- exp(fits$theta)
#         fm_theta <- surveillance::plapply(seq_len(S), apply_optimise_spp_theta,
#                                           first_fit, fits,
#                                           G, disty, pis,
#                                           .parallel = control$cores,
#                                           .verbose = FALSE)
#         theta <- unlist(lapply(fm_theta, `[[`, 1))
#         theta <- log(theta)
#         fits$theta <- update_coefs(log(fits$theta),theta,control$update_kappa[2])
#         }
#     }
#     # e-step - get the log-likes and taus
#     logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
#     taus <- get_taus(pis, logls_mus$logl_sp, G, S)
#     if(ite <3)taus <- shrink_taus(taus,G=G)
#
#
#     #update the likelihood
#     logl_old <- logl_new
#     logl_new <- get_incomplete_logl_sam(eta = additive_logistic(pis,inv = TRUE)[-G],
#                                         first_fit, fits, spp_weights, G, S, disty)
#
#     if(!control$quiet)message(paste0("Iteration ",ite,"\n"))
#     if(!control$quiet)message(paste0("Loglike: ", logl_new,"\n"))
#     if(!control$quiet)message(c("Pis: ", paste(pis," ")," \n"))
#
#     ite <- ite + 1
#   }
#
#   taus <- data.frame(taus)
#   names(taus) <- paste("grp.", 1:G, sep = "")
#   names(pis) <- paste("G", 1:G, sep = ".")
#   eta <- additive_logistic(pis, TRUE)[-G]
#
#   return(list(logl = logl_new, alpha = fits$alpha, beta = fits$beta,
#               gamma = fits$gamma, delta = fits$delta, theta = fits$theta,
#               eta = eta, pis = pis, taus = round(taus,4),
#               first_fit = first_fit))
#
# }

# "initiate_fit_sam" <- function(y, X, W, U, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, size, control){
#
#
#   # if(control$init_glmnet)
#     fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_glmnet_sam_inits, y, X, W, U,
#                                        site_spp_weights, offset, y_is_na, disty, size,
#                                       .parallel = control$cores, .verbose = FALSE)
#   # else fm_sp_mods <-  surveillance::plapply(seq_len(S), apply_glm_sam_inits, y, X, W, U,
#   #                                                             site_spp_weights, offset, y_is_na, disty, size,
#   #                                                             .parallel = control$cores, .verbose = FALSE)
#
#
#   alpha <- unlist(lapply(fm_sp_mods, `[[`, 1))
#   if(ncol(X)==1){
#     beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))[,-1,drop=FALSE]
#     } else {
#     beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))
#   }
#   if(ncol(X)==0) beta <- rep(-999999,G)
#   if(ncol(W)>1){
#     gamma <- do.call(rbind,lapply(fm_sp_mods, `[[`, 3))
#   } else {
#     gamma <- unlist(lapply(fm_sp_mods, `[[`, 3))
#   }
#   if(!is.null(U)){
#     if(ncol(U)>1){
#       delta_spp <- do.call(rbind,lapply(fm_sp_mods, `[[`, 4))
#     } else {
#       delta_spp <- matrix(unlist(lapply(fm_sp_mods, `[[`, 4)),ncol=1)
#     }
#   } else {
#     delta_spp <- matrix(unlist(lapply(fm_sp_mods, `[[`, 4)),ncol=1)
#   }
#   theta <- unlist(lapply(fm_sp_mods, `[[`, 5))
#
#   # species_to_remove <- which(apply(beta, 1, function(x) all(is.na(x))))
#   #
#   # if(length(species_to_remove)>0){
#   #   #update fits
#   #   alpha <- alpha[-species_to_remove]
#   #   beta <- beta[-species_to_remove,,drop=FALSE]
#   #   theta <- theta[-species_to_remove]
#   #
#   #   # update y, y_is_na and weights
#   #   # updated_y <- update_species_data_structure(y, y_is_na, spp_weights, site_spp_weights, species_to_remove)
#   #   # y <- updated_y[[1]]
#   #   # y_is_na <- updated_y[[2]]
#   #   # site_spp_weights <- updated_y[[3]]
#   # } else {
#   #   species_to_remove <- NA
#   # }
#
#   prev_min_sites <- control$minimum_sites_occurrence
#   sel_omit_spp <- which(colSums(y>0, na.rm = TRUE) <= prev_min_sites)
#   if(length(sel_omit_spp)>0) beta <- beta[-sel_omit_spp,,drop=FALSE]
#
#   if(G==1) control$init_method <- 'kmeans'
#
#   if(control$init_method=='kmeans'){
#     # if(!control$quiet)message( "Initial groups by K-means clustering\n")
#     fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
#     tmp_grp <- fmmvnorm$cluster
#     grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
#     grp_coefs <- matrix(grp_coefs,nrow=G)
#   }
#
#   if(control$init_method=='kmed'){
#   # message( "Initial groups parameter estimates by K-medoids\n")
#       mrwdist <- kmed::distNumeric(beta, beta, method = "mrw")
#       fmmvnorm <- kmed::fastkmed(mrwdist, ncluster = G, iterate = 100)
#       tmp_grp <- fmmvnorm$cluster
#       grp_coefs <- beta[fmmvnorm$medoid,,drop=FALSE]
#       grp_coefs <- matrix(grp_coefs,nrow=G)
#   }
#
#   if(control$init_method=='random2'){
#     fmmvnorm <- stats::kmeans(beta, centers=G, nstart=100)
#     tmp_grp <- fmmvnorm$cluster
#     grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
#     grp_coefs <- matrix(grp_coefs,nrow=G) # not check this...
#
#     random_coefs <- sam_random_inits(alpha, grp_coefs, gamma, delta_spp, theta, S, G, X, W, U, disty, mult=0.3)
#
#     alpha <- random_coefs[[1]]
#     grp_coefs <- random_coefs[[2]]
#     gamma <- random_coefs[[3]]
#     # delta <- random_coefs[[4]]
#     if(disty%in%c(4,6))theta <- random_coefs[[5]]
#   }
#
#  if(ncol(X[,,drop=FALSE])==1){
#    names(grp_coefs)[2] <- names(X[,,drop=FALSE])[2]
#  } else {
#    colnames(grp_coefs) <- colnames(X[,,drop=FALSE])
#  }
#
#   ## get the means for delta all coeff.
#   delta <- apply(delta_spp,2,mean,trim=0.1) #trim away extreme values
#   names(delta) <- colnames(U)
#
#   #get taus as starting values
#   if(G==1){
#     taus <- matrix(1,nrow=ncol(y), ncol = G)
#   } else {
#   taus <- matrix(0,nrow=ncol(y), ncol= G)
#   if(length(sel_omit_spp)>0){
#     for(j in 1:length((1:S)[-sel_omit_spp]))
#       taus[(1:S)[-sel_omit_spp][j],fmmvnorm$cluster[j]] <- 1
#       taus[sel_omit_spp,] <- matrix(runif(length(sel_omit_spp)*G),length(sel_omit_spp), G)
#       } else {
#   for(j in seq_len(S))taus[j,fmmvnorm$cluster[j]] <- 1
#         }
#   }
#
#   taus <- taus/rowSums(taus)
#   taus <- shrink_taus(taus,G=G)
#   pis <- colMeans(taus)
#
#   results <- list()
#   results$grps <- tmp_grp
#   results$alpha <- alpha
#   results$beta <- grp_coefs
#   results$gamma <- gamma
#   results$delta <- delta
#   results$theta <- theta
#   results$taus <- taus
#   results$pis <- pis
#
#   return(results)
# }

"sam_optimise" <- function(y, X, W, U, offset, spp_weights, site_spp_weights, y_is_na,
                           S, G, disty, size, powers, start_vals, control){

  inits <- c(start_vals$alpha, start_vals$beta, start_vals$eta, start_vals$gamma,
             start_vals$delta, start_vals$theta)
  npx <- as.integer(ncol(X))
  if(ncol(W)>1) npw <- as.integer(ncol(W[,-1,drop=FALSE])) else npw <- as.integer(0)
  if(!is.null(U)) npu <- as.integer(ncol(U)) else npu <- as.integer(0)

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
  if( npw > 0){
    gamma.score <- as.numeric(matrix( NA, nrow=S, ncol=ncol(W)))
  } else {
    gamma.score <- -999999
  }
  if( npu > 0){
    delta.score <- as.numeric(matrix(NA, ncol=ncol(U)))
  } else {
    delta.score <- -999999
  }
  if( npw > 0){
    control$optiPart <- as.integer(1)
    gamma.score <- as.numeric(rep(NA, length(gamma)))
  } else {
    control$optiPart <- as.integer(0)
    gamma.score <- -999999
  }
  if(disty%in%c(4,6)){
    control$optiDisp <- as.integer(1)
    theta.score <- as.numeric(rep(NA, length(theta)))
  }else{
    control$optiDisp <- as.integer(0)
    theta.score <- -999999
  }
  scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,
                                       gamma.score,delta.score,theta.score))))
  control$conv_cpp <- as.integer(0)

  #model quantities
  pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, S))
  loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

  if(control$print_cpp_start_vals)print_starting_values(as.integer(S),
                                                        as.integer(G),
                                                        as.integer(npx),
                                                        as.integer(npw),
                                                        as.integer(n),
                                                        as.double(alpha),
                                                        as.double(beta),
                                                        as.double(eta),
                                                        as.double(gamma),
                                                        as.double(theta))

  #c++ call to optimise the model (needs pretty good starting values)
  tmp <- .Call("species_mix_cpp",
               as.numeric(as.matrix(y)), as.numeric(as.matrix(X)), as.numeric(as.matrix(W[,-1,drop=FALSE])), as.numeric(offset),
               as.numeric(spp_weights), as.numeric(as.matrix(site_spp_weights)), as.integer(as.matrix(!y_is_na)), as.numeric(size),
               # SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_weights, SEXP Rsite_spp_weights, SEXP Ry_not_na, // data
               as.integer(S), as.integer(G), as.integer(npx), as.integer(npw), as.integer(n),
               as.integer(disty),as.integer(control$optiDisp),as.integer(control$optiPart),
               # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
               as.double(alpha), as.double(beta), as.double(eta), as.double(gamma), as.double(theta),
               # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
               alpha.score, beta.score, eta.score, gamma.score, theta.score, as.integer(control$getscores), scores,
               # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
               pis_out, mus, loglikeS, loglikeSG,
               # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
               as.integer(control$maxit_cpp), as.integer(control$trace_cpp), as.integer(control$nreport_cpp),
               as.numeric(control$abstol_cpp), as.numeric(control$reltol_cpp), as.integer(control$conv_cpp), as.integer(control$printparams_cpp),
               # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
               as.integer(control$optimise_cpp), as.integer(control$loglOnly_cpp), as.integer(control$derivOnly_cpp),
               # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
               PACKAGE = "ecomix")

  ret <- tmp
  ret$logl <- ret$logl * -1
  ret$mus <- array(mus, dim=c(n, S, G))
  if(npw>0){
    if(!disty%in%c(4,6))
      ret$coefs <- list(alpha = ret$alpha, beta = matrix(ret$beta,G,npx),
                        eta = ret$eta,gamma = matrix(ret$gamma,S,npw))
    else
      ret$coefs <- list(alpha = ret$alpha, beta = matrix(ret$beta,G,npx),
                        eta = ret$eta,gamma = matrix(ret$gamma,S,npw),
                        theta = ret$theta)
  } else {
    if(!disty%in%c(4,6))
      ret$coefs <- list(alpha = ret$alpha, beta = matrix(ret$beta,G,npx),
                        eta = ret$eta)
    else
      ret$coefs <- list(alpha = ret$alpha, beta = matrix(ret$beta,G,npx),
                        eta = ret$eta, theta = ret$theta)
  }
  ret$names <- list(spp=colnames(y), SAMs=paste("SAM", 1:G, sep=""),
                    Xvars=colnames(X), Wvars=colnames(W[,-1,drop=FALSE]))

  if(!disty%in%c(4,6))
    ret$scores <- list(alpha.scores = alpha.score, beta.scores = beta.score,
                       eta.scores=eta.score,gamma.scores = gamma.score)
  else
    ret$scores <- list(alpha.scores = alpha.score, beta.scores = beta.score,
                       eta.scores=eta.score,gamma.scores = gamma.score,
                       theta.scores=theta.score)

  ret$S <- S; ret$G <- G; ret$npx <- npx; ret$npw <- npw; ret$n <- n;
  ret$start.vals <- inits
  ret$loglikeSG <- matrix(loglikeSG,  nrow = S, ncol = G)  #for residuals
  ret$loglikeS <- loglikeS  #for residuals
  ret$removed_species <- start_vals$first_fit$removed_species
  gc()
  return(ret)
}




###### SAM internal functions ######

"shrink_taus" <- function(taus, G){
  if( G==1)
    return( taus)
  magical.alpha <- (1-0.8*G)/(0.8*(2-G)-1) ## Dunstan et al., 2013, JABES
  taus_star <- (2*magical.alpha*taus-magical.alpha+1)/(2*magical.alpha - magical.alpha*G + G)
  return(taus_star)
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
                              y_is_na, size, archetype_formula, species_formula,
                              control, family,removed_species)  {
  if( titbits==TRUE)
    titbits <- list( Y = y, X = X, W = W, U = U, spp_weights = spp_weights,
                     site_spp_weights = site_spp_weights, offset = offset,
                     y_is_na = y_is_na, size = size, archetype_formula =  archetype_formula,
                     species_formula = species_formula, control = control,
                     family = family, removed_species = removed_species)
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
    if( "archetype_formula" %in% titbits)
      titbits$archetype_formula <- archetype_formula
    if( "species_formula" %in% titbits)
      titbits$species_formula <- species_formula
    if( "control" %in% titbits)
      titbits$control <- control
    if( "family" %in% titbits)
      titbits$family <- family
    if( "removed_species" %in% titbits)
      titbits$removed_species <- removed_species
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

"species_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))

  # stopifnot(all(is.finite(x)))
}

"covariate_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))
  # stopifnot(all(is.finite(x)))
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


"update_coefs" <- function(old, new, kappa=1){
  if(any(is.na( old)))
    tmp <- new
  else
    tmp <- old + kappa*(new-old)
  return( tmp)
}


# "update_species_data_structure" <- function(y, y_is_na, spp_weights, site_spp_weights, species_to_remove){
#   if(is.na(species_to_remove)) return(list(y, y_is_na, spp_weights, site_spp_weights))
#   else return(list(y[,-species_to_remove], y_is_na[,-species_to_remove],
#                    spp_weights[-species_to_remove], site_spp_weights[,-species_to_remove]))
# }


"print_starting_values" <-  function(S, G, npx, npw, n, alpha, beta, gamma,
                                     eta, theta){
  message(S, "species.\n")
  message(G, "groups.\n")
  message(npx, "archetype coefs.\n")
  message(npw, "species coefs.\n")
  message(n,"sites.\n")
  message("starting species intercepts:\n",alpha,"\n")
  message("starting archetype parameters:\n",beta,"\n")
  message("starting archetype membership:\n",additive_logistic(eta),"\n")
  message("starting species parameters:\n",gamma,"\n")
  message("starting species specific dispersion parameters:\n",theta,"\n")
}

# "reltol_fun" <- function(logl_n1, logl_n){
#   return(abs(logl_n1 - logl_n) > (abs(logl_n1 - logl_n) / abs(logl_n)))
# }

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


"sam_random_inits" <- function(alpha, beta, gamma, delta, theta, S, G, X, W, U, disty, mult=0.3, control.sd = control$init_sd){
  if(is.null(control.sd)){
    my.sd <- mult*sd( alpha); if( is.na( my.sd)) my.sd <- 0.1
  } else {
    my.sd <- control.sd
  }
  alpha <- alpha + rnorm(S, sd = my.sd)
  if(is.null(control.sd)){
    my.sd <- mult*sd( beta); if( is.na( my.sd) | my.sd==0) my.sd <- control.sd
  } else {
    my.sd <- control.sd
  }
  beta <- beta + as.numeric(matrix(rnorm(G * ncol(X), mean = 0, sd = my.sd), ncol = ncol(X), nrow = G))
  if( ncol(W)>1){
    if(is.null(control.sd)){
      my.sd <- mult*sd(gamma); if(is.na(my.sd)|my.sd==0) my.sd <- control.sd
    } else {
      my.sd <- control.sd
    }
    gamma <- gamma + as.numeric( matrix(rnorm(S*ncol(W[,-1,drop=FALSE]), mean=0, my.sd), ncol=ncol(W[,-1,drop=FALSE]), nrow=S))
  }
  if( !is.null(U)){
    if(is.null(control.sd)){
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


"sam_internal_pred_groups" <- function(alpha, beta, taus, gamma, G, S, X, W,
                                       offset = NULL, family){

  if (family %in% c("bernoulli","binomial"))
    link.fun <- make.link("logit")
  if (family %in% c("negative.binomial","poisson","ippm"))
    link.fun <- make.link("log")
  if (family %in% "gaussian")
    link.fun <- make.link("identity")
  if (is.null(offset))
    offset <- rep(0, nrow(X))

  outpred_arch <- matrix(NA, dim(X)[1], G)
  colnames(outpred_arch) <- paste("G", 1:G, sep = ".")

  for (g in seq_len(G)) {
    s.outpred <- matrix(NA, dim(X)[1], length(alpha))
    for (s in seq_len(S)) {
      etaMix <- as.numeric(X%*%beta[g, ])
      if(ncol(W)>1) etaSpp <- as.numeric(W%*%c(alpha[s],gamma[s, ]))
      else etaSpp <- alpha[s]
      eta <- etaMix + etaSpp + offset
      s.outpred[, s] <- link.fun$linkinv(eta)
    }

    outpred_arch[, g] <- apply(s.outpred*rep(taus[, g],each = dim(X)[1]),
                               1, sum)/sum(taus[, g])
  }
  return(outpred_arch)
}

"sam_internal_pred_species" <- function(alpha, beta, taus, gamma, G, S, X, W,
                                        offset = NULL, family){

  if (family %in% c("bernoulli","binomial"))
    link.fun <- make.link("logit")
  if (family %in% c("negative.binomial","poisson","ippm"))
    link.fun <- make.link("log")
  if (family %in% "gaussian")
    link.fun <- make.link("identity")
  if (is.null(offset))
    offset <- rep(0, nrow(X))

  outpred_spp <- matrix(0, dim(X)[1], S)

  for (g in seq_len(G)) {
    etaMix <- matrix(as.numeric(X%*%beta[g, ]), nrow(X), S, byrow=FALSE)
    if(ncol(W)>1) etaSpp <- W%*%t(cbind(alpha,gamma))
    else etaSpp <- matrix(alpha, nrow(X), S, byrow=TRUE)
    eta <- etaMix + etaSpp + offset
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

  if(is.list(inits)){
    alpha <- as.numeric(inits$alpha)
    beta <- as.numeric(inits$beta)
    if(!is.null(npw))
      gamma <- as.numeric(inits$gamma)
    else
      gamma <- -99999
    eta <- as.numeric(inits$eta)
    if(npw>0){
      gamma <- as.numeric(inits$gamma)
    } else {
      gamma <- -999999
    }
    if(disty%in%c(4,6)){
      theta <- as.numeric(inits$theta)
    } else {
      theta <- -999999
    }
    if(return_list) res <- list(alpha=alpha,beta=beta,eta=eta,gamma=gamma,theta=theta)
    else res <- c(alpha,beta,eta,gamma,theta)
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
      gamma <- -999999
      # start <- start + S
    }
    if(disty%in%c(4,6)){
      theta <- inits[start + 1:S]
    } else {
      theta <- -999999
    }

    if(return_list) res <- list(alpha=alpha,beta=beta,eta=eta,gamma=gamma,theta=theta)
    else res <- c(alpha,beta,eta,gamma,theta)
  }

  return(res)

}

"set_control_sam" <- function(control){
  if(is.null(control))control <- species_mix.control()
  control
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

"calc_info_crit_sam" <-  function(tmp) {

    # K <- length(tmp$beta) + tmp$G + tmp$S + ifelse(ncol(W)>1,ncol(W)-1,0) + ifelse(!is.null(U),ncol(U),0) + ifelse(disty %in% c(4,5,6),S,0)

    k <- length(unlist(tmp$coefs))
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


"plapply" <- function (X, FUN, ..., .parallel = 1, .seed = NULL, .verbose = TRUE) {
  if (!(useCluster <- inherits(.parallel, "cluster"))) {
    stopifnot(length(.parallel) == 1L, is.vector(.parallel,
                                                 "numeric"), .parallel >= 1)
    .parallel <- as.vector(.parallel, mode = "integer")
    if (.Platform$OS.type == "windows" && .parallel > 1L) {
      useCluster <- TRUE
      .parallel <- parallel::makeCluster(.parallel)
      on.exit(parallel::stopCluster(.parallel))
    }
  }
  FUN <- match.fun(FUN)
  .FUN <- if (useCluster || is.primitive(FUN)) {
    FUN
  }
  else {
    verboseExpr <- if (isTRUE(.verbose)) {
      if (.parallel == 1L && interactive()) {
        env <- new.env(hash = FALSE, parent = environment(FUN))
        environment(FUN) <- env
        env$pb <- txtProgressBar(min = 0, max = length(X),
                                 initial = 0, style = 3)
        on.exit(close(env$pb), add = TRUE)
        quote(setTxtProgressBar(pb, pb$getVal() + 1L))
      }
      else {
        on.exit(cat("\n"), add = TRUE)
        quote(cat("."))
      }
    }
    else if (is.call(.verbose) || is.expression(.verbose)) {
      .verbose
    }
    else if (is.character(.verbose)) {
      on.exit(cat("\n"), add = TRUE)
      substitute(cat(.verbose))
    }
    do.call(add.on.exit, list(FUN, verboseExpr))
  }
  if (!is.null(.seed)) {
    if (useCluster) {
      parallel::clusterSetRNGStream(cl = .parallel, iseed = .seed)
    }
    else {
      if (!exists(".Random.seed", envir = .GlobalEnv,
                  inherits = FALSE)) {
        set.seed(NULL)
      }
      .orig.seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .orig.seed, envir = .GlobalEnv),
              add = TRUE)
      if (.parallel == 1L) {
        set.seed(seed = .seed)
      }
      else {
        stopifnot(requireNamespace("parallel", quietly = TRUE))
        set.seed(seed = .seed, kind = "L'Ecuyer-CMRG")
        parallel::mc.reset.stream()
      }
    }
  }
  if (useCluster) {
    parallel::parLapply(cl = .parallel, X = X, fun = .FUN,
                        ...)
  }
  else if (.parallel == 1L) {
    lapply(X = X, FUN = .FUN, ...)
  }
  else {
    parallel::mclapply(X = X, FUN = .FUN, ..., mc.preschedule = TRUE,
                       mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = .parallel)
  }
}

"add.on.exit" <- function (FUN, expr){
  FUN <- match.fun(FUN)
  if (is.null(expr <- substitute(expr))) {
    return(FUN)
  }
  if (is.primitive(FUN)) {
    stop("not implemented for primitive functions")
  }
  onexitexpr <- substitute(on.exit(expr))
  obody <- body(FUN)
  body(FUN) <- if (is.call(obody) && identical(as.name("{"),
                                               obody[[1L]])) {
    as.call(append(x = as.list(obody), values = onexitexpr,
                   after = 1L))
  }
  else {
    as.call(c(as.name("{"), onexitexpr, obody))
  }
  FUN
}


