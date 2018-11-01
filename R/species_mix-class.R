##### Species mix functions to export ####

#' @title species_mix objects
#' @rdname species_mix
#' @name species_mix
#' @description Fits a finite mixture model to identify species archetype
#' models (SAMs).
#' @details species_mix is used to fit mixtures of glms to multivariate
#' species data. The function uses BFGS to optimise the mixture likelihood.
#' There is the option to use EM algorithm to get appropriate starting parameters.
#' `species_mix` acts as a wrapper for fitmix.cpp that allows for easier data
#' input. The data frames are merged into the appropriate format for the use
#' in fitmix.cpp. Minima is found using vmmin (BFGS). Currently 'bernoulli',
#' 'poisson', 'ippm' (inhomogenous Poisson point process), 'negative_binomial'
#'  and 'tweedie' distributions can be fitted using the species_mix function.
#' @param archetype_formula an object of class "formula" (or an object that can be
#' coerced to that class). The response variable (left hand side of the
#' formula) needs to be either 'presence', 'occurrence', 'abundance',
#' 'biomass' or 'quantity' data. The type of reponse data will help specify
#' the type of error distribution to be used. The dependent variables
#' (the right hind side) of this formula specifies the dependence of the
#' species archetype probabilities on covariates. For all model the basic
#' formula structure follows something like this:
#' cbind(spp1,spp2,spp3)~1+temperature+rainfall
#' @param species_formula an object of class "formula" (or an object that can be
#' coerced to that class). The right hand side of this formula specifies the
#' dependence of the species"'" data on covariates (typically different covariates
#' to \code{archetype_formula} to avoid confusing confounding). Current the formula
#' is set at ~ 1 by default for species-specific intercepts for the archetype models.
#' If you include a species specific formula which has more than an intercept you
#' will be fitting a partial species archetype model which has species
#' specific covariates and archetype specific covariates.
#' @param data a matrix of dataframe which contains the 'species_data'
#' matrix, a const and the covariates in the strucute of spp1, spp2, spp3,
#' const, temperature, rainfall. dims of matirx should be
#' nsites*(nspecies+const+covariates).
#' @param n_mixtures The number of mixing components (groups) to fit.
#' @param distribution The family of statistical distribution to use within
#' the ecomix models. a  choice between "bernoulli", "poisson",
#' "ippm" (inhomogeneous Poisson point process model), "negative_binomial", "tweedie"
#' and "gaussian" distributions are possible and applicable to specific types
#' of data.
#' @param offset a numeric vector of length nrow(data) (n sites) that is included into
#' the model as an offset. It is included into the conditional part of the model
#' where conditioning is performed on the SAM.
#' @param weights a numeric vector of length ncol(Y) (n species) that is used as weights
#' in the log-likelihood calculations. If NULL (default) then all weights are
#' assumed to be identically 1. Because we are estimating the log-likelihood
#' over species (rather than sites), the weights should be a vector n species
#' long. The exception is under the use of the 'ippm' distribution where
#' weights must be a nrow(data)*n_species matrix, which provides a
#' species-specific background weights used to estimate the species-specific
#' marginal likelihoods.
#' @param bb_weights a numeric vector of n species long. This is used for undertaking
#' a Bayesian Bootstrap. See 'vcov.species_mix' for more details.
#' @param control a list of control parameters for optimisation and calculation.
#' See details. From \code{species_mix.control} for details on optimistaion
#' parameters.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficents. These are distribution specific, but at a
#' minimum you will need pis (additive_logitic transformed), alpha
#' (intercepts) and beta (mixing coefs).
#' @param standardise Booliean. If TRUE, standarise the covariate data.
#' @param titbits either a boolean or a vector of characters. If TRUE (default for species_mix(qv)), then some objects used in the estimation of the model"'"s parameters are returned in a list entitled "titbits" in the model object. Some functions, for example plot.regimix(qv) and predict.regimix(qv), will require some or all of these pieces of information. If titbits=FALSE (default for species_mix.multifit(qv)), then an empty list is returned. If a character vector, then just those objects are returned. Possible values are:"Y" for the outcome matrix, "X" for the model matrix for the RCP model, "offset" for the offset in the model, "site_spp_weights" for the model weights, "archetype_formula" for the formula for the SAMs, "species_formula" for the formula for the species-specific model, "control" for the control arguments used in model fitting, "dist" for the conditional distribution of the species data. Care needs to be taken when using titbits=TRUE in species_mix.multifit(qv) calls as titbits is created for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then setting titbits=TRUE may give users problems with memory.
#' @importFrom graphics abline hist legend lines matplot par plot points polygon rect
#' @importFrom stats as.formula binomial cooks.distance cov cutree dbinom dist dnbinom dnorm dpois
#' fitted gaussian glm hclust lm logLik model.matrix model.offset model.response
#' model.weights pbinom pnbinom pnorm poisson
#' ppois predict qnorm qqnorm quantile rbinom
#' residuals rgamma rnbinom rnorm rpois runif
#' sd uniroot update update.formula
#' @export
#' @examples
#' library(ecomix)
#' set.seed(42)
#' sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))
#' sp_form <- ~ 1
#' theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,1.9),3,3,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
#' dat[,-1] <- scale(dat[,-1])
#' simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
#'                                             dat,theta,dist="bernoulli")
#' data <- make_mixture_data(species_data = simulated_data$species_data,
#'                                 covariate_data = simulated_data$covariate_data[,-1])
#' fm1 <- species_mix(sam_form, sp_form, data, distribution = 'bernoulli',
#'  n_mixtures=3)

"species_mix" <- function(archetype_formula = NULL, species_formula = stats::as.formula(~1), data,
                          n_mixtures = 3, distribution="bernoulli", offset=NULL,
                          weights=NULL, bb_weights=NULL, control=NULL, inits=NULL,
                          standardise = TRUE, titbits = TRUE){

  data <- as.data.frame(data)
  #the control parameters
  control <- set_control_sam(control)
  if(!control$quiet)
    message( "SAM modelling")
  call <- match.call()
  if(!is.null(archetype_formula))
    archetype_formula <- stats::as.formula(archetype_formula)
  else{
    if(!control$quiet)
      message("There is no SAM model! Please provide a model (intercept at least) -- exitting now")
    return(NULL)
  }
  if(!is.null(species_formula))
    species_formula <- stats::as.formula(species_formula)

  # Create model matrix
  mf <- match.call(expand.dots = FALSE)
  if(distribution=="ippm"){
    m <- match(c("data","offset"), names(mf), 0L)
  } else {
    m <- match(c("data","offset","weights"), names(mf), 0L)
  }
  # m <- match(c("data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  if(distribution=="ippm"){
    mf$na.action <- "na.pass"
  } else {
    mf$na.action <- "na.exclude"
  }
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))

  # get the model matrix and find the fitting formula.
  dat <- clean_data_sam(mf, archetype_formula, NULL, distribution)

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
    message( "There are ", n_mixtures, " archtypes to group the species into")

  # get archetype model matrix
  X <- get_X_sam(archetype_formula, dat$mf.X)

  # get species model matrix
  # W <- get_W_sam(species_formula, dat$mf.W) # don't need yet. but will be important for partial sams.

  #get distribution
  disty.cases <- c("bernoulli","poisson","ippm","negative_binomial","tweedie","gaussian")
  disty <- get_distribution_sam(disty.cases, distribution)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names,distribution)
  spp_weights <- check_spp_weights(bb_weights,S)

  # cat(colnames(site_spp_weights),"\n")

  if(distribution=='ippm'){
    if(!all(colnames(y)==colnames(site_spp_weights))){
      stop(cat('When modelling a inhomogenous poisson point process model,\n species data colnames must match weights colnames.\n\nSpecies data colnames from "data" are:\n',colnames(y),'.\n\nWhile the colnames of the weights are:\n', colnames(site_spp_weights),'\n'))
    }
    if(any(dim(y)!=dim(site_spp_weights))){
      stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
    }
  }

  s.means <- NULL
  s.sds <- NULL
  if (standardise == TRUE) {
    stand.X <- standardise.X(X[, -1])
    X <- as.matrix(cbind(1, stand.X$X))
    s.means <- stand.X$dat.means
    s.sds <- stand.X$dat.sds
  }

  # summarising data to console
  print_input_sam(y, X, S, archetype_formula, species_formula, distribution, quiet=control$quiet)

  # fit this bad boy. bad boys, bad boys, what you gonna do when they come for you.
  tmp <- species_mix.fit(y=y, X=X, G=n_mixtures, S=S, spp_weights=spp_weights,
                         site_spp_weights=site_spp_weights,
                         offset=offset, disty=disty, y_is_na=y_is_na,
                         control=control, inits=inits)

  tmp$dist <- disty.cases[disty]


  tmp$pis <- additive_logistic(tmp$eta)

  #calc posterior porbs and pis.
  if(n_mixtures>1)
    tmp$taus <- calc_post_probs_sam(tmp$pis,tmp$loglikeSG)

  tmp$pis <- colSums(tmp$taus)/S

  #Information criteria
  tmp <- calc_info_crit_sam(tmp)

  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_sam(titbits, y, X, spp_weights, site_spp_weights, offset,
                                 y_is_na , archetype_formula, species_formula,
                                 control, disty.cases[disty])
  class(tmp) <- c("species_mix")
  return(tmp)
}

# @rdname species_mix
# @name species_mix.fit
# @param y is a matrix genertated from \link[stats]{model.response} containing the species information. The matrix has the dimensions n_sites * n_species.
# @param X is a design matrix for the archetype_formula dimension n_sites * n_covariates.
# @param W is a design matrix for species_formula and will be implemented if species_formula has covariates.
# @param G is the number of species archetypes that are being estimated.
# @param S is the number of species to be modelled (this will be calculated internally in species_mix())
# @param spp_weights These are weights on the species logls and are specifically used in the Bayesian Boostrap.
# @param site_spp_weights These are site and species specific weights. For most distributions these will be the same across all species. But this form is required to correctly estiamte the IPPMs. See \link[ecomix]{species_mix} for more details.
# @param offset this is a vector of site specific offsets, this might be something like area sampled at sites.
# @param y_is_na This is a logical matrix used specifically with 'ippm' modelling - don't worry about this, it'll be worked out for you. Yay!
# @param disty the error distribution to used in species_mix estimation. Currently, 'bernoulli', 'poisson', 'ippm' (Poisson point process), 'negative_binomial' and 'guassian' are avaliable - internal conversion of distribution to a integer.
# @param control this is a list of control parameters that alter the specifics of model fitting. See \link[ecomix]{species_mix.control} for details.
# @param inits This will be a vector of starting values for species_mix (i.e you've fitted a model and want to refit it).
# @export

"species_mix.fit" <- function(y, X, G, S, spp_weights, site_spp_weights, offset, y_is_na=NULL, disty, control, inits=NULL){

  disty.cases <- c("bernoulli","poisson","ippm",
                   "negative_binomial","tweedie","gaussian")
  distribution <- disty.cases[disty]
  if(!any(distribution==c("bernoulli","poisson","ippm",
                          "negative_binomial","tweedie","gaussian")))
    stop('current only please check the distribution you are fitting')
  # sp.form <- update(archetype_formula,obs~1+.)

  # need to insert two new calls:
  # get_starting_values_sam() which will generate starting values based on EM or clustering of coefs.
  # sam_optimise() which will fit the model in cpp
  # S <- ncol(y) # how many species are there?

  if(is.null(inits)){

    starting_values  <-  ecomix:::get_starting_values_sam(y = y, X = X,
                                                 spp_weights = spp_weights,
                                                 site_spp_weights = site_spp_weights,
                                                 offset = offset,
                                                 y_is_na = y_is_na,
                                                 G = G, S = S,
                                                 disty = disty,
                                                 control = control)

  } else {
    if(!control$quiet)message('Be careful! You are using your own initial starting values to optimise the species_mix model.')
    inits <- setup_inits_sam(inits, S=S, G=G, np=ncol(X[,-1]), disty, return_list = TRUE)
    print(inits)
    starting_values <- inits
  }

  tmp <- sam_optimise(y=starting_values$first_fit$y,
                      X=starting_values$first_fit$x,
                      offset = starting_values$first_fit$offset,
                      spp_weights = starting_values$first_fit$spp_weights,
                      site_spp_weights = starting_values$first_fit$site_spp_weights,
                      y_is_na = starting_values$first_fit$y_is_na,
                      S = ncol(starting_values$first_fit$y),
                      G = G,
                      Obs = nrow(y),
                      disty = disty,
                      start_vals = starting_values,
                      control = control)

  return(tmp)
}

#'@rdname species_mix
#'@name species_mix.fit
#'@param nstart for species_mix.multifit only. The number of random starts to perform for re-fitting. Default is 10, which will need increasing for serious use.
#'@param mc.cores for species_mix.multifit only. The number of cores to spread the re-fitting over.
#'@export
#'@examples
#' \dontrun{
#' fmods <- species_mix.multifit(sam_form, sp_form, data, distribution = 'bernoulli', nstart = 10, n_mixtures=3)
#' }
"species_mix.multifit" <- function(archetype_formula = NULL, species_formula = stats::as.formula(~1), data,
           n_mixtures = 3, nstart = 10, mc.cores=1, distribution="bernoulli", offset=NULL,
           weights=NULL, bb_weights=NULL, control=species_mix.control(), inits=NULL,
           standardise = TRUE, titbits = TRUE){

    #the control parameters
    control <- set_control_sam(control)
    if(!control$quiet)
      message( "SAM modelling")
    call <- match.call()
    if(!is.null(archetype_formula))
      archetype_formula <- stats::as.formula(archetype_formula)
    else{
      if(!control$quiet)
        message("There is no SAM model! Please provide a model (intercept at least) -- exitting now")
      return(NULL)
    }
    if(!is.null(species_formula))
      species_formula <- stats::as.formula(species_formula)

    # Create model matrix
    mf <- match.call(expand.dots = FALSE)
    if(distribution=="ippm"){
      m <- match(c("data","offset"), names(mf), 0L)
    } else {
      m <- match(c("data","offset","weights"), names(mf), 0L)
    }
    # m <- match(c("data","offset","weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    if(distribution=="ippm"){
      mf$na.action <- "na.pass"
    } else {
      mf$na.action <- "na.exclude"
    }
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    # need this for the na.omit step
    rownames(mf)<-seq_len(nrow(mf))

    # clean the data
    dat <- clean_data_sam(mf, archetype_formula, NULL, distribution)

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
      message( "There are ", n_mixtures, " archtypes to group the species into")

    # get archetype model matrix
    X <- get_X_sam(archetype_formula, dat$mf.X)

    # get species model matrix
    # W <- get_W_sam(species_formula, dat$mf.W) # don't need yet. but will be important for partial sams.

    #get distribution
    disty.cases <- c("bernoulli","poisson","ippm","negative_binomial","tweedie","gaussian")
    disty <- get_distribution_sam(disty.cases, distribution)

    # get offsets
    offset <- get_offset_sam(dat$mf.X)

    # get the weights
    species_names <- colnames(y)
    site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names,distribution)
    spp_weights <- check_spp_weights(bb_weights,S)

    # cat(colnames(site_spp_weights),"\n")

    if(distribution=='ippm'){
      if(!all(colnames(y)==colnames(site_spp_weights))){
        cat(colnames(y),"\n")
        cat(colnames(site_spp_weights),"\n")
        stop(cat('When modelling a inhomogenous poisson point process model,\n species data colnames must match weights colnames.\n\nSpecies data colnames from "data" are:\n',colnames(y),'.\n\nWhile the colnames of the weights are:\n', colnames(site_spp_weights),'\n'))
      }
      if(any(dim(y)!=dim(site_spp_weights))){
        stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
      }
    }

    # standardise the data
    s.means <- NULL
    s.sds <- NULL
    if (standardise == TRUE) {
      stand.X <- standardise.X(X[, -1])
      X <- as.matrix(cbind(1, stand.X$X))
      s.means <- stand.X$dat.means
      s.sds <- stand.X$dat.sds
    }

    # summarising data to console
    print_input_sam(y, X, S, archetype_formula, species_formula, distribution, quiet=control$quiet)

  tmp_fun <- function(x){
      if( !control$quiet & nstart>1)
        setTxtProgressBar(pb, x)
      tmpQuiet <- control$quiet
      control$quiet <- TRUE
      tmp <- species_mix.fit(y=y, X=X, G=n_mixtures, S=S, spp_weights=spp_weights,
                             site_spp_weights=site_spp_weights,
                             offset=offset, disty=disty, y_is_na=y_is_na,
                             control=control, inits=inits)

      tmp$dist <- disty.cases[disty]


      tmp$pis <- additive_logistic(tmp$eta)

      #calc posterior porbs and pis.
      if(n_mixtures>1)
        tmp$taus <- calc_post_probs_sam(tmp$pis,tmp$loglikeSG)

      tmp$pis <- colSums(tmp$taus)/S

      #Information criteria
      tmp <- calc_info_crit_sam(tmp)

      #titbits object, if wanted/needed.
      tmp$titbits <- get_titbits_sam(titbits, y, X, spp_weights, site_spp_weights, offset,
                                     y_is_na , archetype_formula, species_formula,
                                     control, disty.cases[disty])
      class(tmp) <- c("species_mix", distribution)
      return( tmp)
   }

  #    require( parallel)
  if( !control$quiet & nstart>1)
    pb <- txtProgressBar(min = 1, max = nstart, style = 3, char = "-^^,--,~ ")

   #Fit the model many times
   many_starts <- surveillance::plapply(seq_len(nstart), tmp_fun, .parallel = mc.cores, .verbose = !control$quiet)
   return(many_starts)
}


#'@rdname species_mix
#'@name control
#'@param quite Should any reporting be performed? Default is FALSE, for reporting.
#'@param trace int 1=model will report parameter estimates and loglikelihood at each iteration. 0=quite.
#'@param reltol function that determines the relative tolernace for model convergence. Default is quite strict.
#'@param maxit Maximum number of evaluations of the objective function allowed. Defaults to 500.
#'@param cores The number of cores to use in fitting of species mix models. These will be largely used to model the species-specific parameteres.
#'@param em_prefit Logical if TRUE the model will run a slower EM algorithim fit to find starting values.
#'@param em_steps int Default is 3, the number of EM iterations to get to starting values.
#'@param em_refit int Default is 1, number of times to refit using EM.

#'@export
"species_mix.control" <- function(maxit = 1000,
                                  quiet = FALSE,
                                  trace = 1,
                                  cores = 1,
                                  ## intialisation controls
                                  init_method = 'kmeans',
                                  init_sd = 1,
                                  minimum_sites_prevelance = 0.05,
                                  # minimum_occurrence_tolerance_ippm = 20,
                                  ## EM algorithim controls
                                  em_prefit = TRUE,
                                  em_steps = 3,
                                  em_refit = 1,
                                  # em_maxit = 3,
                                  em_abstol = sqrt(.Machine$double.eps),
                                  em_reltol = reltol_fun,
                                  em_maxtau = 0.8,
                                  ## c++ controls
                                  maxit_cpp = 1000,
                                  trace_cpp = 1,
                                  nreport_cpp = 1,
                                  abstol_cpp = sqrt(.Machine$double.eps),
                                  reltol_cpp = sqrt(.Machine$double.eps),
                                  conv_cpp = 1,
                                  printparams_cpp = 0,
                                  optimise_cpp = 1,
                                  loglOnly_cpp = 0,
                                  derivOnly_cpp = 0,
                                  getscores_cpp = 0, ...){
               #general controls
  rval <- list(maxit = maxit, quiet = quiet, trace = trace,
               cores = cores,
               #initialisation controls
               init_method = init_method, init_sd = init_sd, minimum_sites_prevelance = minimum_sites_prevelance,
               # minimum_occurrence_tolerance_ippm = minimum_occurrence_tolerance_ippm,
               #em controls
               em_prefit = em_prefit, em_refit = em_refit, em_steps = em_steps,# em_maxit = em_maxit,
               em_abstol = em_abstol, em_reltol = em_reltol, em_maxtau = em_maxtau,
               #cpp controls
               maxit_cpp = maxit_cpp, trace_cpp = trace_cpp, nreport_cpp = nreport_cpp,
               abstol_cpp = abstol_cpp, reltol_cpp = reltol_cpp, conv_cpp = conv_cpp,
               printparams_cpp = printparams_cpp, optimise_cpp = optimise_cpp,
               loglOnly_cpp = loglOnly_cpp, derivOnly_cpp = derivOnly_cpp,
               getscores_cpp = getscores_cpp)
  rval <- c(rval, list(...))
  if (is.null(rval$em_reltol))
    rval$em_reltol <- sqrt(.Machine$double.eps)
  rval
}

#' @rdname species_mix
#' @name simulate_species_mix_data
#' @param archtype_formula formula to simulate species_mix data, needs to have the format: cbind(spp1,spp2,spp3,...,sppN)~1 + x1 + x2
#' @param species_formula formula to simulate species_mix species-specific responses, e.g: ~1
#' @param dat a matrix of variables to simulate data from.
#' @param theta coefficents for each species archetype. Matrix of G x number of parameters. Each row is a different species archetype.
#' @param distribution Which statistical distribution to simulate data for. 'bernoulli', 'gaussian', 'ippm', 'negative_binomial','poisson' and 'tweedie'.
#' @export
#' @examples
#' \dontrun{
#' archetype_formula <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))
#' theta <- matrix(c(-0.9,-0.6,0.5,1,-0.9,1,0.9,-0.9,2.9,-1,0.2,-0.4),4,3,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
#' simulated_data <- simulate_species_mix_data(archetype_formula,~1,dat,theta,dist="bernoulli")
#' }
## need to update this to take the new formula framework and simulate ippm data.
"simulate_species_mix_data" <-  function (archetype_formula, species_formula, dat, theta, distribution = "bernoulli"){

  if(distribution=='ippm'){#stop('simulation of ippm data has not been set up yet, watch this space')
  message('simulating Inhomogenous Point Process Data on a regular 100x100 grid')

    n_sp <- length(archetype_formula[[2]])-1
    n_g <- dim(theta)[1]
    x <- y <- 1:100 / 100
    grid2D <- expand.grid( x, y)
    grid2D$cellArea <- rep( 1/100, nrow( grid2D))  #all cells have same size here
    grid2D$x1 <- runif(nrow(grid2D))
    grid2D$x2 <- runif(rnorm(grid2D))

    sp_name <- all.vars(archetype_formula)[seq_len(n_sp)]

    X <- as.matrix(data.frame(const=1,x1=grid2D$x1,x2=grid2D$x2))
    lambdas <- matrix(0, dim(X)[1], n_sp, dimnames=list(NULL,sp_name))
    sp_int <- rep(0, n_sp)
    group <- rep(0, n_sp)
    for (s in seq_len(n_sp)) {
      g <- sample(n_g,1)
      sp_int[s] <- runif(1, -1, .5)
      log_lambda <-  X%*%c(sp_int[s],theta[g,-1])
      lambdas[, s] <- exp(log_lambda)
      group[s] <- g
    }

    LAMBDAS <- apply(lambdas,2,function(x)sum(x*grid2D$cellArea))
    Ns <- sapply(LAMBDAS,function(x)rpois(n=1, lambda= x))  #the observed number of presences
    preds_df <- data.frame(idx=1:nrow(X),X)
    presences <- list()
    for(i in seq_len(n_sp)){
      presences[[i]] <- sample(x=preds_df$idx,size=Ns[i], replace=TRUE, prob=lambdas[,i]/LAMBDAS[i])# TRUE
    }

    presence_coords <- lapply(presences,function(x)grid2D[x,1:2])
    presences_sort <- lapply(presences,sort)

    sp_dat_po_ul<- data.frame(sp=rep(sp_name,unlist(lapply(presences,length))),cell_num=unlist(presences_sort))
    po_matrix <- table_to_species_data(sp_dat_po_ul,site_id = 'cell_num',species_id = 'sp')
    po_matrix[po_matrix==0]<-NA
    po_covariates <- X[as.numeric(rownames(po_matrix)),]
    presence_data <- data.frame(po_matrix,po_covariates)
    bkdata <- cbind(matrix(0,nrow(X),n_sp),X)
    colnames(bkdata) <- colnames(presence_data)
    mm <- rbind(presence_data,bkdata)
    dat <- mm[c(sp_name,"const","x1","x2")]

    species_specific_cell_counts <- lapply(seq_along(sp_name),
                                           function(x)table(sp_dat_po_ul[sp_dat_po_ul$sp==sp_name[x],2]))

    df <- data.frame(id=preds_df$idx,area=grid2D$cellArea,x1=grid2D$x1)

    sp_weights <- lapply(seq_along(sp_name),
                         function(x)(weights=df$area/as.numeric(species_specific_cell_counts[[x]][match(df$id,as.numeric(names(species_specific_cell_counts[[x]])))])))

    sp_weights_mat <- data.frame(cell_id = 1:10000, do.call(cbind,sp_weights))

    m <- sp_weights_mat
    presence_sites <- m[rowSums(is.na(m[,-1]))!=ncol(m[,-1]), ]
    presence_sites <- data.frame(presence_sites)
    background_sites <- data.frame(cell_id=1:10000,matrix(rep(grid2D$cellArea,n_sp),
                                                          nrow(grid2D),n_sp))

    wts <- rbind(presence_sites[,-1],background_sites[,-1])
    colnames(wts) <- c(sp_name)#,"const","x1","x2")
    y <- dat[,1:n_sp]
    X <- dat[,c(n_sp+1):ncol(dat)]
    y_is_na <- is.na(y)
    weights <- wts
    offset <- rep(0,nrow(dat))
    pi <- tapply(group, group, length)/n_sp
    return(list(species_data = y, covariate_data = X, background_weights = wts, offset=offset,
         y_is_na=y_is_na,group = group, pi = pi, sp.int = sp_int, lambdas = LAMBDAS))

  } else {

  S <- length(archetype_formula[[2]])-1
  #update the formula to old format.
  if(!is.null(archetype_formula))
    archetype_formula <- stats::as.formula(archetype_formula)
  form_org <- archetype_formula
  form <- update(archetype_formula,y~.)

  #check the species formula.
  sp_form <- species_formula
  species_int_coefs <- check_species_formula(sp_form)

  X <- stats::model.matrix(form, dat)
  out <- matrix(0, dim(X)[1], S)
  k <- dim(theta)[1]
  if(dim(theta)[2]!=ncol(X))stop('theta must have the same dimensions as "data" (do not forget the intercept)')
  sp.int <- rep(0, S)
  group <- rep(0, S)
  for (s in 1:S) {
    g <- ceiling(stats::runif(1) * k)
    if (distribution == "bernoulli") {
      theta[g, 1] <- stats::runif(1, -3, 3)
      sp.int[s] <- theta[g, 1]
      lgtp <- X %*% theta[g, ]
      p <- exp(lgtp)/(1 + exp(lgtp))
      out[, s] <- rbinom(dim(X)[1], 1, p)
    }
    if (distribution == "negative_binomial") {
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = TRUE) > 5000 | sum(tmp) < 100) {
        theta[g, 1] <- stats::runif(1, -15, 5)
        sp.int[s] <- theta[g, 1]
        lgtp <- X %*% theta[g, ]
        p <- exp(lgtp)
        tmp <- rnbinom(dim(X)[1], mu = p, size = 1)
      }
      out[, s] <- tmp
    }
    if (distribution == "poisson") {
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = TRUE) > 5000 | sum(tmp) < 100) {
        theta[g, 1] <- stats::runif(1, -5, 5)
        sp.int[s] <- theta[g, 1]
        lgtp <- X %*% theta[g, ]
        tmp <- rpois(dim(X)[1], lambda = exp(lgtp))
      }
      out[, s] <- tmp
    }
    # if (distribution == "tweedie") {
    #   tmp <- rep(6e+05, dim(X)[1])
    #   while (max(tmp, na.rm = TRUE) > 5e+05 | sum(tmp) < 100) {
    #     theta[g, 1] <- stats::runif(1, -15, 5)
    #     theta[g, 1] <- stats::runif(1, 1, 5)
    #     sp.int[s] <- theta[g, 1]
    #     lgtp <- X %*% theta[g, ]
    #     p <- exp(lgtp)
    #     tmp <- rTweedie(dim(X)[1], mu = p, phi = 2, p = 1.6)
    #   }
    #   out[, s] <- tmp
    # }
    if (distribution == "gaussian") {
      sp.int <- NULL
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = TRUE) > 50000 | sum(tmp) < 100) {
        theta[g, 1] <- stats::runif(1, 1, 500)
        sp.int[s] <- theta[g, 1]
        # theta[g, 1] <- stats::runif(1, 1, 500)
        lgtp <- X %*% theta[g, ]
        p <- (lgtp)
        tmp <- stats::rnorm(dim(X)[1], mean = p, sd = 1)
      }
      out[, s] <- tmp
    }
    group[s] <- g
  }
  pi <- tapply(group, group, length)/S
  colnames(out) <- all.vars(form_org)[1:S]
  out <- as.matrix(out)
  dat <- as.matrix(dat)
  return(list(species_data = out, covariate_data = dat, group = group, pi = pi, sp.int = sp.int))
  }
}

#'@rdname species_mix
#'
#'@name species_mix_estimate_groups
#'
#'@description This function runs the 'species_mix' function but it iterates through groups (1 to 10) by default.
#'
#'@export
#'
#'@examples
#'
#' \dontrun{
#' fm_groups <- species_mix_estimate_groups(sam_form, sp_form, data, distribution = 'bernoulli', n_mixtures=1:10)
#'}
"species_mix_estimate_groups" <- function(archetype_formula = NULL, species_formula = ~1,
                                          data, n_mixtures = 1:10, distribution="bernoulli",
                                          offset=NULL, weights=NULL, control= species_mix.control(cores = 1),
                                          inits=NULL, standardise = TRUE, titbits=FALSE){
  my.fun <- function(G, archetype_formula, species_formula, data, distribution, offset,
                     weights, bb_weights, control, inits, standardise, titbits){

    message("Fitting group",G,"\n")
    tmp <- species_mix(archetype_formula, species_formula, data, G, distribution, offset,
                        weights, bb_weights, control, inits, standardise, titbits)
    return(tmp)
    }

  out <- surveillance::plapply(n_mixtures, my.fun, archetype_formula, species_formula, data, distribution, offset,
                               weights, bb_weights, control, inits, standardise, titbits,
                               .parallel = control$cores, .verbose = !control$quiet)
  aic <- rep(0,length(n_mixtures))
  bic <- rep(0,length(n_mixtures))
  fm <- list()
  for(i in seq_along(n_mixtures))
    if(!is.atomic(out[[i]])){
      aic[i] <- out[[i]]$aic
      bic[i] <- out[[i]]$bic
      fm[[i]] <- list(logl=out[[i]]$logl,coef=out[[i]]$coef,tau=out[[i]]$taus,pi=out[[i]]$pi)#,covar=out[[i]]$covar)
    }
  return(list(aic=aic,bic=bic,fm=fm))
}

#' @rdname species_mix
#'
#' @export
#'
"coef.species_mix" <- function (object, ...){
  res <- list()
  res$alpha <- object$coefs$alpha
  # names( res$alpha) <- object$names$spp
  if( !is.null( object$coef$beta)){
    res$beta <- matrix(object$coefs$beta, nrow = object$G, ncol = object$np)
    # colnames( res$tau) <- object$names$spp # need to keep the spp names
  }
  if(!is.null( object$coef$disp)){
    res$logDisp <- object$coef$disp
    # names( res$logDisp) <- object$names$spp
  }
  return(res)
}

#'@rdname species_mix
#'
#'@export
#'
#'@examples
#'
#'#Print information about a species_mix model
#'\dontrun{
#'print(fm1)}

"print.species_mix" <-  function (x,...){
  cat(x$titbits$distribution, "species_mix model\n")
  cat("\nMixing probabilities\n")
  print(x$pi)
  cat("\nCoefficents\n")
  print(x$coef)
  # if(!is.na(x$se[1])){
  #   cat("\nStandard Errors of coefficents\n")
  #   print(x$se)
  # }
  cat("\nPosterior Probabilities\n")
  print(x$taus)
}

#' @rdname species_mix
#' @export
#' @description  The randomised quantile residuals ("RQR", from Dunn and Smyth, 1996) are defined by their marginal distribution function (marginality is over #' other species observations within that site; see Foster et al, in prep).

"residuals.species_mix" <- function( object, ..., type="RQR", quiet=FALSE) {
    if( ! type %in% c("deviance","RQR"))
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
      switch( object$dist,
              bernoulli = { fn <- function(y,mu,logdisp,power) pbinom( q=y, size=1, prob=mu, lower.tail=TRUE)},
              poisson = { fn <- function(y,mu,logdisp,power) ppois( q=y, lambda=mu, lower.tail=TRUE)},
              ippm = { fn <- function(y,mu,logdisp,power) ppois( q=y, lambda=mu, lower.tail=TRUE)},
              negative_binomial = { fn <- function(y,mu,logdisp,power) pnbinom( q=y, mu=mu, size=1/exp( logdisp), lower.tail=TRUE)},
              gaussian = { fn <- function(y,mu,logdisp,power) pnorm( q=y, mean=mu, sd=exp( logdisp), lower.tail=TRUE)})

      for( ss in 1:object$S){
        if( all( object$titbits$power==-999999))  tmpPow <- NULL else tmpPow <- object$titbits$power[ss]
        if( object$dist %in% c("bernoulli","poisson","negative_binomial")){
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
        if( object$dist == "gaussian"){
          tmp <- fn( object$titbits$Y[,ss], object$mus[,ss,], object$coef$disp[ss], object$titbits$power[ss])
          tmp <- rowSums( tmp * object$pis)
          resids[,ss] <- qnorm( tmp)
        }
      }
      if( !quiet & sum( resids==Inf | resids==-Inf)>0)
        message( "Some residuals, well",sum( resids==Inf | resids==-Inf), "to be precise, are very large (infinite actually).\nThese observations lie right on the edge of the realistic range of the model for the data (maybe even over the edge).")

    }
    if( type=="RQR.sim"){
      nsim <- 1000
      if( is.null( mc.cores))
        mc.cores <- getOption("mc.cores", 4)
      resids <- matrix( NA, nrow=object$n, ncol=object$S)
      RQR.fun <- function(ii){
        if( !quiet)
          setTxtProgressBar(pb, ii)
        X1 <- kronecker( matrix( 1, ncol=1, nrow=nsim), fm$titbits$X[ii,,drop=FALSE])
        W1 <- kronecker( matrix( 1, ncol=1, nrow=nsim), fm$titbits$W[ii,,drop=FALSE])
        sims <- simRCPdata( nRCP=object$nRCP, S=object$S, n=nsim, p.x=object$p.x, p.w=object$p.w, alpha=object$coef$alpha, tau=object$coef$tau, beta=object$coef$beta, gamma=object$coef$gamma, logDisps=object$coef$disp, powers=object$titbits$power, X=X1, W=W1, offset=object$titbits$offset,dist=object$dist)
        sims <- sims[,1:object$S]
        yi <- object$titbits$Y[ii,,drop=FALSE]
        many_yi <- matrix( rep( yi, each=nsim), ncol=object$S)
        F_i <- colMeans( sims <= many_yi)
        F_i_minus <- colMeans( sims < many_yi)
        r_i <- runif( object$S, min=F_i_minus, max=F_i)
        return( qnorm( r_i))
      }
      if( !quiet)
        pb <- txtProgressBar(min = 1, max = object$n, style = 3, char = "><(('> ")
      if( Sys.info()['sysname'] == "Windows" | mc.cores==1)
        resids <- lapply( 1:object$n, RQR.fun)
      else
        resids <- parallel::mclapply( 1:object$n, RQR.fun, mc.cores=mc.cores)
      if( !quiet)
        message("")
      resids <- matrix( unlist( resids), nrow=object$n, ncol=object$S, byrow=TRUE)
      if( !quiet & sum( resids==Inf | resids==-Inf)>0)
        message( "Some residuals, well",sum( resids==Inf | resids==-Inf), "to be precise, are very large (infinite actually).\nThese observations lie right on the edge of the Monte Carlo approximation to the distribution function.\nThis may be remedied by getting a better approximation (increasing nsim).")
    }
    return( resids)
  }


#'@rdname species_mix
#'
#'@export
#'
#'@examples
#'
#'# Estimate the variance-covariance matrix.
#'# This will provide estimates of uncertainty for model parameters.
#'\dontrun{
#' vcov(fm1)}
"vcov.species_mix" <- function (object, ..., object2=NULL, method = "FiniteDifference",
                                nboot = 10, mc.cores = 1, D.accuracy=2){
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
    p.x <- ncol(X[,-1])
    offy <- object$titbits$offset
    spp_wts <- object$titbits$spp_weights
    site_spp_wts <- object$titbits$site_spp_weights
    Y <- object$titbits$Y
    y_is_na <- object$titbits$y_is_na
    distribution <- object$titbits$distribution
    disty.cases <- c("bernoulli","poisson","ippm","negative_binomial","tweedie","gaussian")
    disty <- get_distribution_sam(disty.cases, distribution)
    S <- object$S
    G <- object$G
    n <- object$n
    control <- object$titbits$control

    # values for optimisation.
    inits <- object$coefs
    np <- as.integer(ncol(X[,-1]))
    n <- Obs <- as.integer(nrow(X))
    start_vals <- setup_inits_sam(inits,S,G,np,disty,return_list = TRUE)

    # parameters to optimise
    alpha <- as.numeric(start_vals$alpha)
    beta <- as.numeric(start_vals$beta)
    eta <- as.numeric(start_vals$eta)
    disp <- as.numeric(start_vals$disp)

    #scores
    alpha.score <- as.numeric(rep(NA, length(alpha)))
    beta.score <- as.numeric(rep(NA, length(beta)))
    eta.score <- as.numeric(rep(NA, length(eta)))
    disp.score <- as.numeric(rep(NA, length(disp)))
    getscores <- 1
    scores <- as.numeric(rep(NA,length(c(alpha,beta,eta,disp))))


    if(disty%in%c(4,6)){
      control$optiDisp <- as.integer(1)
    }else{
      control$optiDisp <- as.integer(0)
    }

    #model quantities
    pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
    mus <- as.numeric(array( NA, dim=c( Obs, S, G)))  #container for the fitted spp model
    loglikeS <- as.numeric(rep(NA, S))
    loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

    #c++ call to optimise the model (needs pretty good starting values)

    if (method %in% c("FiniteDifference")) {
      grad_fun <- function(x) {
        #x is a vector of first order derivates to optimise using numDeriv in order to find second order derivates.
        start <- 0
        alpha <- x[start + seq_len(S)]
        start <- start + S
        beta <- x[start + seq_len((G*np))]
        start <- start + (G*np)
        eta <- x[start + seq_len(G - 1)]
        start <- start + (G-1)
        if(disty%in%c(4,6))
          disp <- x[start + seq_len(S)]
        else
          disp <- rep(-999999,S)
        tmp <- .Call("species_mix_cpp",
                 as.numeric(as.matrix(Y)), as.numeric(as.matrix(X[,-1])), as.numeric(offy), as.numeric(spp_wts),
                 as.numeric(as.matrix(site_spp_wts)), as.integer(as.matrix(!y_is_na)),
                 # SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_weights, SEXP Rsite_spp_weights, SEXP Ry_not_na, // data
                 as.integer(S), as.integer(G), as.integer(np), as.integer(Obs), as.integer(disty),
                 as.integer(control$optiDisp),
                 # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
                 as.double(alpha), as.double(beta), as.double(eta), as.double(disp),
                 # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
                 alpha.score, beta.score, eta.score, disp.score, as.integer(control$getscores), scores,
                 # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
                 pis_out, mus, loglikeS, loglikeSG,
                 # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
                 as.integer(control$maxit_cpp), as.integer(control$trace_cpp), as.integer(control$nreport_cpp),
                 as.numeric(control$abstol_cpp), as.numeric(control$reltol_cpp), as.integer(control$conv_cpp),
                 as.integer(control$printparams_cpp),
                 # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
                 as.integer(0), as.integer(0), as.integer(1),
                 # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly,
                 PACKAGE = "ecomix")

        tmp1 <- c(alpha.score, beta.score, eta.score)
        # if( p.w > 0)#class( object$titbits$species_formula) == "formula")
          # tmp1 <- c( tmp1, gamma.score)
        if(any(!is.na(object$coef$disp)|!is.null(object$coef$disp)))
          tmp1 <- c( tmp1, disp.score)
        return(tmp1)
      }
      mod_coefs <- unlist(object$coefs)
      # if(any(is.na( object$coef$disp))) mod_coefs <- mod_coefs[!is.na(mod_coefs)]

      hess <- numDeriv::jacobian(grad_fun, mod_coefs)
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
        coefMat <- sam_bootstrap(object, nboot=nboot, type=method, mc.cores=mc.cores, quiet=TRUE, orderSamps=FALSE)
      else
        coefMat <- object2
      vcov.mat <- cov( coefMat)
    }
    return(vcov.mat)
  }

#' @rdname species_mix
#' @export
"AIC.species_mix" <- function (object, ..., k = 2){
  p <- length(unlist(object$coefs))
  if (is.null(k))
    k <- 2
  star.ic <- -2 * object$logl + k * p
  return(star.ic)
}

#' @rdname species_mix
#' @export
"BIC.species_mix" <-  function (object, ...){
    p <- length(unlist(object$coefs))
    k <- log(object$n)
    star.ic <- -2 * object$logl + k * p
    return(star.ic)
  }

#' @rdname species_mix
#' @export
"summary.species_mix" <-function (object, ...){
    if (is.null(object$vcov)) {
      object$vcov <- matrix(NA, nrow = length(unlist(object$coef)),
                            ncol = length(unlist(object$coef)))
      stop("No variance matrix has been supplied")

    }
    res <- cbind(unlist(object$coefs), sqrt(diag(object$vcov)))
    res <- cbind(res, res[, 1]/res[, 2])
    res <- cbind(res, 2 * (1 - pnorm(abs(res[, 3]))))
    colnames(res) <- c("Estimate", "SE", "z-score", "p")
    return(res)
  }

#'@rdname species_mix
#'@param object is a matrix model returned from the species_mix model.
#'@param newobs a matrix of new observations for prediction.
#'@export
#'@examples
#'\dontrun{
#'fm1 <- species_mix(form,data)
#'preds_fm1 <- predict(fm1,newdata)}

"predict.species_mix" <- function(object, newobs=NULL, ...) {

   family <- object$dist

  if(is.null(newobs)) newobs <- object$titbits$X[,-1]
  if(ncol(newobs) != ncol(object$coefs$beta)) stop("Number of coefficients does not match the number of predictors in the new observations - double check you're not including an intercept.")

  G <- object$G
  S <- object$S
  n <- object$n

  ## To calculate fitted values, take a linear combination of the fitted mus from each component
  predict_mus <- predict_y <- predict_etas <- matrix(0,nrow(newobs),S)
  for(gg in seq_len(G)) {
    etas.K <- matrix(object$coefs$alpha,nrow(newobs),S,byrow=TRUE) + matrix(as.matrix(newobs)%*%object$coef$beta[gg,],nrow(newobs),S,byrow=FALSE)
    if(family == "gaussian") mus.K <- etas.K
    if(family == "bernoulli") mus.K <- exp(etas.K)/(1+exp(etas.K))
    if(family %in% c("poisson","ippm","negative_binomial")) mus.K <- exp(etas.K)
    predict_mus <- predict_mus + matrix(object$taus[,gg],nrow(newobs),S,byrow=T)*mus.K
  }

  ## To generate an actual response matrix; simulate component-label then simulated response conditional
  for(ss in seq_len(S)) {
    get.z <- sample(seq_len(G), 1, prob = object$taus[ss,])
    etas.K <- object$alpha[ss] + as.matrix(newobs)%*%object$coefs$beta[get.z,]
    if(family == "bernoulli")
      predict_y[,ss] <- rbinom(n, 1, p=exp(etas.K)/(1+exp(etas.K)))
    if(family == "poisson")
      predict_y[,ss] <- rpois(n, lambda = exp(etas.K))
    if(family == "negative.binomial")
      predict_y[,ss] <- rnbinom(n, mu = exp(etas.K), size = exp(-object$disp[ss]))
  }

  rownames(predict_mus) <- rownames(predict_y) <- 1:nrow(newobs);
  colnames(predict_mus) <- colnames(predict_y) <- names(object$alpha)

  return(list(predict_mus = round(predict_mus,5), predict_y = predict_y))
}

# "predict.species_mix" <-function (object, new_obs, ...){
#   mixture.model <- object
#   if (class(mixture.model)[2] == "bernoulli") {
#     G <- length(mixture.model$pi)
#     covar <- mixture.model$covar[-(1:(G - 1)), -(1:(G -
#         1))]
#     coef <- mixture.model$coef
#     model.fm <- stats::as.formula(mixture.model$formula)
#     model.fm[[2]] <- NULL
#     X <- stats::model.matrix(model.fm, new_obs)
#     link.fun <- stats::make.link("logit")
#     outvar <- matrix(NA, dim(X)[1], G)
#     outpred <- matrix(NA, dim(X)[1], G)
#     colnames(outvar) <- colnames(outpred) <- paste("G", 1:G, sep = ".")
#     for (g in 1:G) {
#       lp <- as.numeric(X %*% coef[g, ])
#       outpred[, g] <- link.fun$linkinv(lp)
#       dhdB <- (exp(lp)/(1 + exp(lp))) * X - exp(lp)^2/((1 +
#           exp(lp))^2) * X
#       c2 <- covar[seq(g, dim(covar)[1], G), seq(g, dim(covar)[1],
#         G)]
#       for (k in 1:dim(X)[1]) {
#         outvar[k, g] <- (dhdB[k, ] %*% c2) %*% (dhdB[k, ])
#       }
#     }
#   }
#   if (class(mixture.model)[2] == "negative_binomial" | class(mixture.model)[2] ==
#       "tweedie") {
#     G <- length(mixture.model$pi)
#     covar <- mixture.model$covar[-1 * c(1:(G - 1), (dim(mixture.model$covar)[1] -
#         length(mixture.model$sp.intercept) + 1):dim(mixture.model$covar)[1]),
#       -1 * c(1:(G - 1), (dim(mixture.model$covar)[1] -
#           length(mixture.model$sp.intercept) + 1):dim(mixture.model$covar)[1])]
#     sp.int <- mixture.model$sp.intercept
#     coef <- mixture.model$coef
#     model.fm <- stats::as.formula(mixture.model$formula)
#     model.fm[[2]] <- NULL
#     X <- cbind(stats::model.matrix(model.fm, new_obs), 1)
#     offset <- stats::model.frame(model.fm, data = new_obs)
#     offset <- stats::model.offset(offset)
#     if (is.null(offset))
#       offset <- rep(0, nrow(X))
#     outvar <- matrix(NA, dim(X)[1], G)
#     outpred <- matrix(NA, dim(X)[1], G)
#     colnames(outvar) <- colnames(outpred) <- paste("G",
#       1:G, sep = ".")
#     for (g in 1:G) {
#       s.outvar <- matrix(NA, dim(X)[1], length(sp.int))
#       s.outpred <- matrix(NA, dim(X)[1], length(sp.int))
#       for (s in seq_along(sp.int)) {
#         lp <- as.numeric(X %*% c(coef[g, ], sp.int[s]) +
#             offset)
#         s.outpred[, s] <- exp(lp)
#         dhdB <- exp(lp) * X
#         c2 <- covar[c(seq(g, G * (dim(X)[2] - 1), G),
#           G * (dim(X)[2] - 1) + s), c(seq(g, G * (dim(X)[2] -
#               1), G), G * (dim(X)[2] - 1) + s)]
#         for (k in 1:dim(X)[1]) {
#           s.outvar[k, s] <- (dhdB[k, ] %*% c2) %*% (dhdB[k,
#             ])
#         }
#       }
#       outpred[, g] <- apply(s.outpred * rep(mixture.model$tau[,
#         g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[,
#           g])
#       outvar[, g] <- apply(s.outvar * rep(mixture.model$tau[,
#         g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[,
#           g])
#     }
#   }
#   if (class(mixture.model)[2] == "ippm" | class(mixture.model)[2] == "poisson") {
#     G <- length(mixture.model$pi)
#     covar <- mixture.model$covar[-1 * c(1:(G - 1)), -1 * c(1:(G - 1))]
#     sp.int <- mixture.model$sp_intercept
#     coef <- mixture.model$coef
#     model.fm <- stats::as.formula(mixture.model$formula)
#     model.fm[[2]] <- NULL
#     X <- cbind(stats::model.matrix(model.fm, new_obs))
#     offset <- stats::model.frame(model.fm, data = new_obs)
#     offset <- stats::model.offset(offset)
#     if (is.null(offset))
#       offset <- rep(0, nrow(X))
#     outvar <- matrix(NA, dim(X)[1], G)
#     outpred <- matrix(NA, dim(X)[1], G)
#     colnames(outvar) <- colnames(outpred) <- paste("G", 1:G, sep = ".")
#     for (g in 1:G) {
#       s.outvar <- matrix(NA, dim(X)[1], length(sp.int))
#       s.outpred <- matrix(NA, dim(X)[1], length(sp.int))
#       for (s in seq_along(sp.int)) {
#         lp <- as.numeric(X %*% c(sp.int[s],coef[g, ]) + offset)
#         s.outpred[, s] <- exp(lp)
#         dhdB <- exp(lp) * X
#         ## pull out the sp intercept and mixing component covariates.
#         c2 <- covar[c(G * (dim(X)[2] - 1) + s,seq(g, G * (dim(X)[2] - 1), G)), c((G * (dim(X)[2] - 1) + s),seq(g, G * (dim(X)[2] - 1), G))]
#         for (k in 1:dim(X)[1]) {
#           s.outvar[k, s] <- (dhdB[k, ] %*% c2) %*% (dhdB[k, ])
#         }
#       }
#       outpred[, g] <- apply(s.outpred * rep(mixture.model$tau[, g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[, g])
#       outvar[, g] <- apply(s.outvar * rep(mixture.model$tau[, g], each = dim(X)[1]), 1, mean)/sum(mixture.model$tau[, g])
#     }
#   }
#   list(fit = outpred, se.fit = sqrt(outvar))
# }

###### SAM internal functions for fitting ######
# replace this function with one that include eta (linear predictor) as an offset in when estimating the species-specific intercepts.
## update the dispersion parameters if needed.
## my own little function to predict from a glm.fit object.
## glmfit object
## newmatrix is a set of new observations or the original obs.
## offset is an offset
## disty is the distribution

"predict.glm.fit" <- function(glmfit, newmatrix, offset, disty){

  if(disty == 1)
    fam <- binomial()
  if(disty == 2 | disty == 3 | disty == 4)
    fam <- poisson()
  if(disty == 6)
    fam <- gaussian()


  coefs <- as.matrix(glmfit$coef)
  eta <- as.numeric(as.matrix(newmatrix) %*% as.numeric(coefs)) + offset
  preds <- fam$linkinv(eta)
  return(preds)
}

## this will give the species intercepts with respect to the mixture linear predictor.
"apply_glm_sam_sp_params" <- function(ss, y, X, G, site_spp_weights, offset,
                                      y_is_na, disty, fits){

  if(disty == 1)
    fam <- binomial()
  if(disty == 2 | disty == 3 | disty == 4)
    fam <- poisson()
  if(disty == 6)
    fam <- gaussian()


  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.numeric(y[ids_i,ss])
  }
  out1 <- kronecker(rep( 1, G), outcomes)
  X1 <- kronecker(rep( 1, G), X[ids_i,])
  wts1 <- kronecker(rep( 1, G), as.numeric(site_spp_weights[ids_i,ss]))
  offy1 <- kronecker(rep( 1, G), offset[ids_i])
  offy2 <- X[ids_i,-1] %*% t(fits$beta)
  offy2 <- as.numeric(offy2)
  offy <- offy1 + offy2

  if( disty != 5){ #don't use for tweedie
  ft_sp <- suppressWarnings(stats::glm.fit(x=as.data.frame(X1),
                            y=as.numeric(out1),
                            weights=as.numeric(wts1),
                            offset=as.numeric(offy),
                            family=fam))
    my_coefs <- coef(ft_sp)
  }
  disp <- NA
  if( disty == 4){
    preds <- predict.glm.fit(ft_sp, X1, offy, disty)
    tmp <- MASS::theta.mm(out1, preds,
                          weights=c(wts1),
                          dfr=length(out1),
                          eps=1e-4)
    if(tmp>2) tmp <- 2
    disp <- log(1/tmp)
  }
  if( disty == 6){
    preds <- predict.glm.fit(ft_sp, X1, offy, disty)
    disp <- log(sqrt(sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard Deviation.
  }

  return(list(alpha = my_coefs[1], beta = my_coefs[-1], disp = disp))

}

## function for starting values.
"apply_glm_sam_inits" <- function(ss, y, X, site_spp_weights, offset, y_is_na, disty){

  # which family to use?
  if(disty == 1)
    fam <- binomial()
  if(disty == 2 | disty == 3 | disty == 4)
    fam <- poisson()
  if(disty == 6)
    fam <- gaussian()

  ids_i <- !y_is_na[,ss]

  if (disty==3){
    outcomes <- as.matrix(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  } else {
    outcomes <- as.matrix(y[ids_i,ss])
  }

  if( disty != 5){
    ft_sp <- try(suppressWarnings(stats::glm.fit(x=as.data.frame(X[ids_i,]),
                            y=as.numeric(outcomes),
                            weights=as.numeric(site_spp_weights[ids_i,ss]),
                            offset=offset[ids_i],
                            family=fam)), silent=TRUE)
    if (class(ft_sp) %in% 'try-error'){
      my_coefs <- rep(NA, ncol(X[ids_i,]))
    } else {
      my_coefs <- coef(ft_sp)
    }
  }
  disp <- NA
  if(disty == 4){
    preds <- predict.glm.fit(ft_sp, X[ids_i,], offset[ids_i], disty)
    tmp <- MASS::theta.mm(outcomes, preds,
                          weights=c(site_spp_weights[ids_i,ss]),
                          dfr=length(outcomes),
                          eps=1e-4)
    if(tmp>2) tmp <- 2
    disp <- log(1/tmp)
  }
  if( disty == 6){
    preds <- predict.glm.fit(ft_sp, X1, offy, disty)
    disp <- log(sqrt(sum((outcomes - preds)^2)/length(outcomes)))  #should be something like the resid standard Deviation.
  }
   return(list(alpha = my_coefs[1], beta = my_coefs[-1], disp = disp))
}


# do I need to include the species intercepts in the offsets?

"apply_glm_group_tau_sam" <- function (gg, y, X, site_spp_weights, offset, y_is_na, disty,
                                       tau, fits, mus){

  ### setup the data stucture for this model.
  Y_tau <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
  X_no_NA <- list()
  for (jj in 1:ncol(y)){
    X_no_NA[[jj]] <- X[!y_is_na[,jj],-1]
  }
  X_tau <- do.call(rbind, X_no_NA)
  n_ys <- sapply(X_no_NA,nrow)

  wts_tau <- rep(tau[,gg],c(n_ys))
  if(disty==4)wts_tau <- rep(tau[,gg],c(n_ys))/(1+rep(exp(-fits$disp),each=n_ys)*as.vector(mus[gg,,]))

  site_weights <- as.matrix(as.matrix(unlist(as.data.frame(site_spp_weights[!y_is_na]))))
  wts_tauXsite_weights <- wts_tau*site_weights
  offy_mat <- replicate(ncol(y),offset)
  offy1 <- unlist(as.data.frame(offy_mat[!y_is_na]))
  offy2 <- fits$alpha[rep(1:length(fits$alpha),n_ys)]
  offy <- as.numeric(offy1 + offy2)

  options(warn = -1)
  # which family to use?
  if( disty == 1)
    fam <- binomial()
  if( disty == 2 | disty == 3 |disty == 4)
    fam <- poisson()
  if( disty == 6)
    fam <- gaussian()
  if (disty==3){
    Y_tau <- as.matrix(Y_tau/site_weights)
  } else {
    Y_tau <- as.matrix(Y_tau)
  }

  if(disty!=5){ #don't use for tweedie
    ft_mix <- suppressWarnings(glm.fit(x = as.data.frame(X_tau),
                      y = as.numeric(Y_tau),
                      weights = c(wts_tauXsite_weights),
                      offset = offy,
                      family = fam))
  }
    mix_coefs <- coef(ft_mix)
    return(as.matrix(mix_coefs))
}

"get_starting_values_sam" <- function(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  temp.warn <- getOption( "warn")
  options( warn=-1)

  #if emfit is in the control do an EM fit to get good starting values for c++
  if(disty==3)control$em_prefit<-FALSE

  if(isTRUE(control$em_prefit)){

    if(!control$quiet)message('Using EM algorithm to find starting values; using',
                              control$em_refit,'refits\n')

    emfits <- list()
    for(ii in seq_len(control$em_refit)){
      emfits[[ii]] <-  fitmix_EM_sam(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control)
    }
    bf <- which.max(vapply(emfits,function(x)c(x$logl),c(logl=0)))
    emfit <- emfits[[bf]]
    start_vals <- list(alpha=emfit$alpha,
                       beta=emfit$beta,
                       disp=emfit$disp,
                       pis=emfit$pis,
                       first_fit = starting_values$first_fit)
  } else {
    if(!control$quiet)message('You are not using the EM algorith to find starting values; starting values are
                              generated using',control$init_method,'\n')
    starting_values <- ecomix:::get_initial_values_sam(y = y, X = X,
                                              spp_weights = spp_weights,
                                              site_spp_weights = site_spp_weights,
                                              offset = offset, y_is_na = y_is_na,
                                              G = G, S = S,
                                              disty=disty,
                                              control = control)
    start_vals <- list(alpha=starting_values$fits$alpha,
                       beta=starting_values$fits$beta,
                       disp=starting_values$fits$disp,
                       pis=starting_values$pis,
                       first_fit = starting_values$first_fit)
  }

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pis, inv = TRUE)[-G]
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)

  return(start_vals)
}

"get_incomplete_logl_sam" <- function(eta, first_fit, fits, spp_weights, G, S, disty){

  pis <- additive_logistic(eta)
  if(is.null(spp_weights))spp_weights <- rep(1,S) #for bayesian boostrap.

  #bernoulli
  if(disty==1){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "logit")
    for(ss in 1:S){
      for(gg in 1:G){
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(lp),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #poisson
  if(disty==2){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(lp),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #ippm
  if(disty==3){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[sp_idx,1] * fits$alpha[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        logl_sp[ss,gg] <- first_fit$y[sp_idx,ss] %*% lp - first_fit$site_spp_weights[sp_idx,ss] %*% exp(lp)
      }
    }
  }
  #negative binomial
  if(disty==4){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss],mu=exp(lp),size=exp(-fits$disp[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  #tweedie
  if(disty==5){
    stop('no tweedie')

  }
  # gaussian
  if(disty==6){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=lp,sd=exp(fits$disp[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  ak <- logl_sp + matrix(rep(log( pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  logl <- sum( sppLogls)

  return(logl)
}

"get_initial_values_sam" <- function(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){

  # get intial model fits
  starting_values <- ecomix:::initiate_fit_sam(y, X, site_spp_weights, offset, y_is_na, G, S, disty, control)

  #if any are errors then remove them from the models for ever.
  updated_y <- ecomix:::update_species_data_structure(y, y_is_na, spp_weights,
                                                      site_spp_weights, starting_values$species_to_remove)
  y <- updated_y[[1]]
  y_is_na <- updated_y[[2]]
  spp_weights <- updated_y[[3]]
  site_spp_weights <- updated_y[[4]]
  S <- ncol(y)

  fits <- list(alpha=starting_values$alpha,beta=starting_values$beta,disp=starting_values$disp)
  first_fit <- list(x = X, y = y, spp_weights = spp_weights,
                    site_spp_weights = site_spp_weights, offset = offset,
                    y_is_na = y_is_na, removed_species = starting_values$species_to_remove)
  logls_mus <- ecomix:::get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
  pis <- rep(1/G, G)
  taus <- ecomix:::get_taus(pis, logls_mus$logl_sp, G, S)
  taus <- ecomix:::skrink_taus(taus, max_tau = 1/G + 0.1, G) #max_tau=1/G + 0.1

  # #now fit the mix model once
  # fmix_coefs <- surveillance::plapply(seq_len(G), ecomix:::apply_glm_group_tau_sam,
  #                                     first_fit$y,
  #                                     first_fit$x,
  #                                     first_fit$site_spp_weights,
  #                                     first_fit$offset,
  #                                     first_fit$y_is_na,
  #                                     disty,
  #                                     taus,
  #                                     fits,
  #                                     logls_mus$fitted,
  #                                     .parallel = control$cores,
  #                                     .verbose = FALSE)#!control$quiet)
  #
  # #update the mix coefs.
  # fmix_coefs <- t(do.call(cbind,fmix_coefs))
  # fits$beta <- ecomix:::update_mix_coefs(fits$beta,fmix_coefs)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$pis <- colMeans(taus)
  res$taus <- taus
  return(res)
}



"get_logls_sam" <- function(first_fit, fits, spp_weights, G, S, disty, get_fitted=TRUE){

  if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(first_fit$y),S))
  # if(is.null(spp_weights))spp_weights <- rep(1,S) #for bayesian boostrap.

  #bernoulli
  if(disty==1){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    link <- stats::make.link(link = "logit")
    for(ss in 1:S){
      for(gg in 1:G){
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        if(get_fitted) fitted_values[gg,,ss] <- link$linkinv(lp)
        logl_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, link$linkinv(lp),log = TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #poisson
  if(disty==2){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        if(get_fitted) fitted_values[gg,,ss] <- exp(lp)
        logl_sp[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(lp),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #ippm
  if(disty==3){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      sp_idx<-!first_fit$y_is_na[,ss]
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        eta <- first_fit$x[sp_idx,1] * fits$alpha[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        # mu <- eta_to_mu(eta,poisson())
        lp <- first_fit$x[sp_idx,1] * fits$alpha[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$beta[gg,] + first_fit$offset[sp_idx]
        if(get_fitted) fitted_values[gg,sp_idx,ss] <- exp(lp)
        # logl_sp[ss,gg] <- sum(first_fit$site_spp_weights[sp_idx,ss] * (first_fit$y[sp_idx,ss] * log(mu) - mu))
        logl_sp[ss,gg] <- (first_fit$y[sp_idx,ss] %*% lp - first_fit$site_spp_weights[sp_idx,ss] %*% exp(lp))
      }
    }
  }

  #negative binomial
  if(disty==4){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        if(get_fitted) fitted_values[gg,,ss] <- exp(lp)
        logl_sp[ss,gg] <- sum(dnbinom(first_fit$y[,ss],mu=exp(lp),size=exp(-fits$disp[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }

  #tweedie
  if(disty==5){
    stop('no tweedie')
  }

  # gaussian
  if(disty==6){
    logl_sp <- matrix(NA, nrow=S, ncol=G)
    for(ss in 1:S){
      for(gg in 1:G){
        #lp is the same as log_lambda (linear predictor)
        lp <- first_fit$x[,1] * fits$alpha[ss] + as.matrix(first_fit$x[,-1]) %*% fits$beta[gg,] + first_fit$offset
        if(get_fitted) fitted_values[gg,,ss] <- lp
        logl_sp[ss,gg] <- sum(dnorm(first_fit$y[,ss],mean=lp,sd=exp(fits$disp[ss]),log=TRUE))
      }
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }
  }
  out.list <- list(logl_sp=logl_sp)
  if(get_fitted) out.list$fitted = fitted_values
  return(out.list)
}

# "eta_to_mu" <- function (eta, family, mu.min = 1e-16, mu.max = 1/mu.min){
#   mu <- family$linkinv(eta)
#   mu[mu < mu.min] = mu.min
#   mu[mu > mu.max] = mu.max
#   mu
# }


"fitmix_EM_sam" <- function(y, X, spp_weights, site_spp_weights, offset, y_is_na, G, S, disty, control){


  ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_sam(y = y, X = X,
                                            spp_weights = spp_weights,
                                            site_spp_weights = site_spp_weights,
                                            offset = offset, y_is_na = y_is_na,
                                            G = G, S = S,
                                            disty=disty,
                                            control = control)

  # first e-step
  fits <- starting_values$fits
  taus <- starting_values$taus
  first_fit <- starting_values$first_fit
  logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){

    # Estimate pis from tau.
    pis <- colSums(taus)/S

    if (any(pis == 0)) {
      starting_values <- get_initial_values_sam(y = y, X = X,
                                                spp_weights = spp_weights,
                                                site_spp_weights = site_spp_weights,
                                                offset = offset, y_is_na = y_is_na,
                                                G = G, S = S,
                                                disty=disty,
                                                control = control)
      pis <- starting_values$pis
      fits <- starting_values$fits
      taus <- starting_values$taus
      first_fit <- starting_values$first_fit
      ite <- 1
    }

    # m-step
    # replace this with a glm.fit version that can deal with the
    # need to include species intercepts and dispersion parameters as an offset in this estimation.
    fmix_coefs <- surveillance::plapply(1:G, apply_glm_group_tau_sam,
                                        y, X, site_spp_weights,
                                        offset, y_is_na, disty, taus,
                                        fits, logls_mus$fitted,
                                        .parallel = control$cores,
                                        .verbose = FALSE)#!control$quiet)

    # update the coefs.
    fmix_coefs_mat <- t(do.call(cbind,fmix_coefs))
    fits$beta <- update_mix_coefs(fits$beta, fmix_coefs_mat)

    fm_sp_int <- surveillance::plapply(1:S, apply_glm_sam_sp_params,
                                       y, X, G, site_spp_weights, offset,
                                       y_is_na, disty, fits,
                                       .parallel = control$cores, .verbose = FALSE)
    #check weights in this.
    alpha <- unlist(lapply(fm_sp_int, `[[`, 1))
    fits$alpha <- update_sp_coefs(fits$alpha,alpha)

    if(disty%in%c(4,6)){
      disp <- unlist(lapply(fm_sp_int, `[[`, 3))
      fits$disp <- update_sp_dispersion(fits$disp,disp)
    }

    # e-step
    # get the log-likes and taus
    logls_mus <- get_logls_sam(first_fit, fits, spp_weights, G, S, disty)
    pis <- rep(1/G, G)
    taus <- get_taus(pis, logls_mus$logl_sp, G, S)

    #update the likelihood
    logl_old <- logl_new
    logl_new <- get_incomplete_logl_sam(eta = additive_logistic(pis,inv = TRUE)[-G], first_fit, fits, spp_weights, G, S, disty)
    ite <- ite + 1#;print(ite)
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", 1:G, sep = "")
  int_out <- fits$alpha
  fm_out <- fits$beta
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- additive_logistic(pis, TRUE)[-1]

  # estimate log-likelihood
  logl_new <- get_incomplete_logl_sam(eta, first_fit, fits, spp_weights, G, S, disty)

  return(list(logl = logl_new, alpha = int_out, beta = fm_out, disp = fits$disp,
              eta = eta, pis = pis, taus = round(taus,4), first_fit = first_fit))

}

"initiate_fit_sam" <- function(y, X, site_spp_weights, offset, y_is_na, G, S, disty, control){

  fm_sp_mods <-  surveillance::plapply(seq_len(S), ecomix:::apply_glm_sam_inits, y, X,
                                       site_spp_weights, offset, y_is_na, disty,
                                      .parallel = control$cores, .verbose = FALSE)

  alpha <- unlist(lapply(fm_sp_mods, `[[`, 1))
  beta <- do.call(rbind,lapply(fm_sp_mods, `[[`, 2))
  disp <- unlist(lapply(fm_sp_mods, `[[`, 3))

  species_to_remove <- which(apply(beta, 1, function(x) all(is.na(x))))
  if(length(species_to_remove)>0){
    #update fits
    alpha <- alpha[-species_to_remove]
    beta <- beta[-species_to_remove,]
    disp <- disp[-species_to_remove]

    # update y, y_is_na and weights
    updated_y <- update_species_data_structure(y, y_is_na, spp_weights, site_spp_weights, species_to_remove)
    y <- updated_y[[1]]
    y_is_na <- updated_y[[2]]
    site_spp_weights <- updated_y[[3]]
  } else {
    species_to_remove <- NA
  }

  if(disty==3){
    n <- max(colSums(y<1,na.rm = TRUE))
    prev_min_sites <- floor(n*control$minimum_sites_prevelance)
    sel_omit_spp <- which(colSums(y>0, na.rm = TRUE) <= prev_min_sites)
  } else {
    n <- nrow(y)
    prev_min_sites <- floor(n*control$minimum_sites_prevelance);
    sel_omit_spp <- which(colSums(y>0, na.rm = TRUE) <= prev_min_sites)
  }

  if(length(sel_omit_spp)>0) beta <- beta[-sel_omit_spp,]

  if(control$init_method=='kmeans'){
    if(disty==3){
      message( "Initial groups parameter estimates by K-medoids\n")
      mrwdist <- kmed::distNumeric(beta, beta, method = "mrw")
      fmmvnorm <- kmed::fastkmed(mrwdist, ncluster = G, iterate = 100)
      tmp_grp <- fmmvnorm$cluster
      grp_coefs <- beta[fmmvnorm$medoid,]
    } else {
    if(!control$quiet)message( "Initial groups by K-means clustering\n")
      tmp1 <- stats::kmeans(beta, centers=G, nstart=100)
      tmp_grp <- tmp1$cluster
      grp_coefs <- apply(beta, 2, function(x) tapply(x, tmp_grp, mean))
    }
  }


  if(control$init_method=='random' | is.null(tmp_grp)){
    if(!control$quiet)message( "Initial groups by random allocation and means from random numbers\n")
    grp_coefs <- matrix( stats::rnorm(G*ncol(beta), sd=control$init_sd, mean=0), nrow=G, ncol=ncol(beta))
    tmp_grp <- sample(1:G, S, replace=TRUE)
  }

  colnames(grp_coefs) <- colnames(X[,-1])
  results <- list()
  results$grps <- tmp_grp
  results$alpha <- alpha
  results$beta <- grp_coefs
  results$disp <- disp
  results$species_to_remove <- species_to_remove

  return(results)
}

"incom_logl_mix_coefs" <- function(x, eta, first_fit, fits, spp_weights, G, S, disty){

  fits$beta <- matrix(x,nrow=nrow(fits$beta),ncol=ncol(fits$beta))
  tmp <- get_incomplete_logl_sam(eta, first_fit, fits, spp_weights, G, S, disty)
  return(-tmp)
}

"sam_optimise" <- function(y, X, offset, spp_weights, site_spp_weights, y_is_na, S, G, Obs, disty, start_vals, control){

  inits <- c(start_vals$alpha, start_vals$beta, start_vals$eta, start_vals$disp)
  np <- as.integer(ncol(X[,-1]))
  n <- as.integer(nrow(X))

  # parameters to optimise
  alpha <- as.numeric(start_vals$alpha)
  beta <- as.numeric(start_vals$beta)
  eta <- as.numeric(start_vals$eta)
  disp <- as.numeric(start_vals$disp)

  #scores
  alpha.score <- as.numeric(rep(NA, length(alpha)))
  beta.score <- as.numeric(rep(NA, length(beta)))
  eta.score <- as.numeric(rep(NA, length(eta)))
  disp.score <- as.numeric(rep(NA, length(disp)))
  getscores <- 1
  scores <- as.numeric(rep(NA,length(c(alpha,beta,eta,disp))))


  if(disty%in%c(4,6)){
    control$optiDisp <- as.integer(1)
  }else{
    control$optiDisp <- as.integer(0)
    # disp <- -99999
  }

  #model quantities
  pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
  mus <- as.numeric(array(NA, dim=c(n, S, G)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, S))
  loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

  #c++ call to optimise the model (needs pretty good starting values)
  tmp <- .Call("species_mix_cpp",
               as.numeric(as.matrix(y)), as.numeric(as.matrix(X[,-1])), as.numeric(offset), as.numeric(spp_weights),
               as.numeric(as.matrix(site_spp_weights)), as.integer(as.matrix(!y_is_na)),
               # SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_weights, SEXP Rsite_spp_weights, SEXP Ry_not_na, // data
               as.integer(S), as.integer(G), as.integer(np), as.integer(n), as.integer(disty),
               as.integer(control$optiDisp),
               # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
               as.double(alpha), as.double(beta), as.double(eta), as.double(disp),
               # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
               alpha.score, beta.score, eta.score, disp.score, as.integer(control$getscores), scores,
               # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
               pis_out, mus, loglikeS, loglikeSG,
               # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
               as.integer(control$maxit_cpp), as.integer(control$trace_cpp), as.integer(control$nreport_cpp),
               as.numeric(control$abstol_cpp), as.numeric(control$reltol_cpp), as.integer(control$conv_cpp),
               as.integer(control$printparams_cpp),
               # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
               as.integer( control$optimise_cpp), as.integer(control$loglOnly_cpp),
               as.integer( control$derivOnly_cpp),
               # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
               PACKAGE = "ecomix")

  ret <- tmp
  ret$logl <- ret$logl * -1
  ret$mus <- array(mus, dim=c(Obs, S, G))

  # beta <- matrix(ret$beta,G,np)
  # colnames(beta) <-

  if(!disty%in%c(4,6))
    ret$coefs <- list(alpha = ret$alpha, beta = matrix(ret$beta,G,np), eta = ret$eta)
  else
    ret$coefs <- list(alpha = ret$alpha, beta = ret$beta, eta = ret$eta, disp = ret$disp)

  ret$names <- list(spp=colnames(y), RCPs=paste("SAM", 1:G, sep=""), Xvars=colnames(X[,-1]))

  if(!disty%in%c(4,6))
    ret$scores <- list(alpha.scores = alpha.score, beta.scores = beta.score, eta.scores=eta.score)
  else
    ret$scores <- list(alpha.scores = alpha.score, beta.scores = beta.score, eta.scores=eta.score, disp.scores=disp.score)

  ret$S <- S; ret$G <- G; ret$np <- np; ret$n <- Obs;
  ret$start.vals <- inits
  ret$loglikeSG <- matrix(loglikeSG,  nrow = S, ncol = G)  #for residuals
  ret$loglikeS <- loglikeS  #for residuals
  return(ret)
}

"sam_bootstrap" <-function (object, nboot=1000, type="BayesBoot", mc.cores=1,
                            quiet=FALSE, orderSamps=FALSE, MLstart=TRUE){
  if (nboot < 1)
    stop( "No Boostrap samples requested.  Please set nboot to something > 1.")
  if( ! type %in% c("BayesBoot","SimpleBoot"))
    stop( "Unknown boostrap type, choices are BayesBoot and SimpleBoot.")
  n.reorder <- 0
  object$titbits$control$optimise <- TRUE #just in case it was turned off
  if(object$titbits$distribution=='ippm')
    stop('IPPM vcov matrix needs to estimated using FiniteDifference method.\n')

  if( type == "SimpleBoot"){
    all.wts <- matrix( sample( 1:object$S, nboot*object$S, replace=TRUE), nrow=nboot, ncol=object$S)
    tmp <- apply( all.wts, 1, table)
    all.wts <- matrix( 0, nrow=nboot, ncol=object$S)
    for( ii in seq_along( tmp))
      all.wts[ii, as.numeric( names( tmp[[ii]]))] <- tmp[[ii]]
  }
  if( type == "BayesBoot")
    all.wts <- object$S * gtools::rdirichlet( nboot, rep( 1, object$S))
  if(MLstart)
    my.inits <- object$coef
  else{
    my.inits <- "random"
    orderSamps <- TRUE
  }

  tmpOldQuiet <- object$titbits$control$quiet
  object$titbits$control$quiet <- TRUE

  my.fun <- function(dummy){
    disty.cases <- c("bernoulli", "poisson", "ippm", "negative_binomial", "tweedie", "gaussian")
    disty <- get_distribution_sam(disty.cases, object$dist)
    # if( !object$titbits$control$quiet){
    # pb <- progress::progress_bar$new(
    #   format = " Running Bayesian bootstrap [:bar] :percent eta: :eta",
    #   total = nboot, clear = FALSE, width= 60)
      # }
    dumbOut <- capture.output(
      samp.object <- species_mix.fit(y=object$titbits$Y,
                                     X=object$titbits$X,
                                     offset = object$titbits$offset,
                                     spp_weights = all.wts[dummy,,drop=TRUE],
                                     site_spp_weights = object$titbits$site_spp_weights,
                                     G = object$G,
                                     S = object$S,
                                     y_is_na = object$titbits$y_is_na,
                                     disty = disty,
                                     control = object$titbits$control,
                                     inits = my.inits))
    # pb$tick()
    if( orderSamps)
      samp.object <- orderPost( samp.object, object)
    return( unlist( samp.object$coef))
  }

  tmp <- surveillance::plapply(seq_len(nboot), my.fun, .parallel = mc.cores)
  boot.estis <- do.call( "rbind", tmp)
  object$titbits$control$quiet <- tmpOldQuiet
  # if( !quiet)
    # message( "")
  # colnames( boot.estis) <- get_long_names_rcp( object)
  class( boot.estis) <- "sam_bootstrap"
  return( boot.estis)
}


###### SAM internal functions ######

"get_taus" <- function(pi, logls, G, S){
  fullLogPis <- matrix(rep(log(pi), each=S), nrow=S, ncol=G)
  a_k <- fullLogPis + logls
  a_m <- apply( a_k, 1, max)
  tmp <- exp( a_k - rep( a_m, times=G))
  log_denom <- a_m + log( rowSums( tmp))
  return( exp( a_k - log_denom))
}

"skrink_taus" <- function( taus, max_tau=0.7, G){
  if( G==1)
    return( taus)
  alpha <- (1-max_tau*G) / ( max_tau*(2-G)-1)
  tau_star <- ( 2*alpha*taus - alpha + 1 ) / ( 2*alpha - alpha*G + G)
  return(tau_star)
}

# "lambda_penalisation_fun" <- function(x,lambda,kappa=0.1){ #assumes that x spans to pretty-well the unpenalised estiamtes
#   min.effective.penalty <- min( which( abs( x-tail( x, 1)) < 0.01 * abs( tail( x, 1))))    #the first that lambda that gives a coef close to the last lambda's corresponding coef
#   min.effective.penalty <- lambda[min.effective.penalty]
#   target.penalty <- kappa * min.effective.penalty
#   res.pos <- which.min( (lambda-target.penalty)^2)
#   res <- x[res.pos]
#   return( res)
# }

"print_input_sam" <- function(y, X, S, archetype_formula, species_formula, distribution, quiet=FALSE){
  if( quiet)
    return( NULL)
  n.tot <- nrow(y)
  if(distribution=='ippm'){
    n_pres <- sum(unlist(y)==1,na.rm=TRUE)
    n_bkgrd <- sum(unlist(y[,1])==0,na.rm=TRUE)
    message("There are ", n_pres, " presence observations for ", S," species")
    message("There are ", n_bkgrd, " background (integration) points for each of the ", S," species")
  } else {
    message("There are ", nrow(X), " site observations for ", S," species")
  }
  archetype_formula[[2]] <- NULL
  message("The model for the SAM is ", Reduce( "paste", deparse(archetype_formula)))
  if(!is.null(species_formula))
  message("The model for the species is ", Reduce( "paste", deparse(species_formula)))
  message("You are implementing a ", distribution, " SAM.")
}

"get_distribution_sam" <- function( disty.cases, dist1) {
  error.msg <- paste( c( "Distribution not implemented. Options are: ", disty.cases, "-- Exitting Now"), collapse=" ")
  disty <- switch( dist1,
                   "bernoulli" = 1, #"bernoulli_sp" = 2, removing bernoulli sp for now as all bernoulli will be species specific ints
                   "poisson" = 2,
                   "ippm" = 3,
                   "negative_binomial" = 4,
                   "tweedie" = 5,
                   "gaussian" = 6,
                   {stop( error.msg)} )
  return( disty)
}

"get_X_sam" <- function(archetype_formula, mf.X){
  form.X <- archetype_formula
  form.X[[2]] <- NULL
  form.X <- stats::as.formula(form.X)
  X <- stats::model.matrix(form.X, mf.X)
  return( X)
}

"get_W_sam" <- function( species_formula, mf.W){
  form.W <- species_formula
  if( !is.null( species_formula)){
    if( length( form.W)>2)
      form.W[[2]] <- NULL #get rid of outcomes
    W <- stats::model.matrix( form.W, mf.W)
    tmp.fun <- function(x){ all( x==1)}
    intercepts <- apply( W, 2, tmp.fun)
    W <- W[,!intercepts,drop=FALSE]
  }
  else
    W <- -999999
  return( W)
}


"species_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))
  # stopifnot(all(is.finite(x)))
}

"covariate_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))
  # stopifnot(all(is.finite(x)))
}

"get_offset_sam"  <- function(mf){
  offset <- stats::model.offset(mf)
  if(any(offset!=0))
    return(offset)
  offset <- rep(0, nrow(mf))
  return(offset)
}

"get_titbits_sam" <- function( titbits, y, X, spp_weights, site_spp_weights, offset,
                               y_is_na , archetype_formula, species_formula,
                               control, distribution)  {
    if( titbits==TRUE)
      titbits <- list( Y = y, X = X, spp_weights=spp_weights,
                       site_spp_weights=site_spp_weights, offset=offset,
                       y_is_na=y_is_na,
                       archetype_formula =  archetype_formula,
                       species_formula = species_formula,
                       control = control,
                       distribution = distribution,
                       removed_species = removed_species)
    else{
      titbits <- list()
      if( "Y" %in% titbits)
        titbits$Y <- y
      if( "X" %in% titbits)
        titbits$X <- X
      # if( "W" %in% titbits)
        # titbits$W <- W
      if( "spp_weights" %in% titbits)
        titbits$spp_weights <- spp_weights
      if( "site_spp_weights" %in% titbits)
        titbits$site_spp_weights <- site_spp_weights
      if( "offset" %in% titbits)
        titbits$offset <- offset
      if( "y_is_na" %in% titbits)
        titbits$y_is_na <- y_is_na
      if( "archetype_formula" %in% titbits)
        titbits$archetype_formula <- archetype_formula
      if( "species_formula" %in% titbits)
        titbits$species_formula <- species_formula
      if( "control" %in% titbits)
        titbits$control <- control
      if( "distribution" %in% titbits)
        titbits$distribution <- distribution
      if( "removed_species" %in% titbits)
        titbits$removed_species <- removed_species
    }
    return( titbits)
}

"get_site_spp_weights_sam"  <- function(mf, site_spp_weights, sp_names, distribution){

  # site_spp_wts <- model.weights(mf)
  if(is.null(site_spp_weights))site_spp_weights <- model.weights(mf)

  if(distribution=='ippm'){
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

"update_mix_coefs" <- function(old, new, kappa=1){
  if(any(is.na( old)))
    tmp <- new
  else
    tmp <- old + kappa*(new-old)
  return( tmp)
}

"update_sp_coefs" <- function(old, new, kappa=1){
  if(any(is.na( old)))
    tmp <- new
  else
    tmp <- old + kappa*(new-old)
  return( tmp)
}

"update_sp_dispersion" <- function(old, new, kappa=1){
  if(any(is.na( old)))
    tmp <- new
  else
    tmp <- old + kappa*(new-old)
  return( tmp)
}

"update_species_data_structure" <- function(y, y_is_na, spp_weights, site_spp_weights, species_to_remove){
  if(is.na(species_to_remove)) return(list(y, y_is_na, spp_weights, site_spp_weights))
  else return(list(y[,-species_to_remove], y_is_na[,-species_to_remove],
                   spp_weights[,-species_to_remove], site_spp_weights[,-species_to_remove]))
}

"reltol_fun" <- function(logl_n1, logl_n){
  return(abs(logl_n1 - logl_n) > (abs(logl_n1 - logl_n) / abs(logl_n)))
}

"setup_inits_sam" <- function(inits, S, G, np, disty, return_list=TRUE){
  if(is.null(inits))res<-NULL

  if(is.list(inits)){
    alpha <- as.numeric(inits$alpha)
    beta <- as.numeric(inits$beta)
    eta <- as.numeric(inits$eta)
    # disp <- as.numeric(inits$disp)
    if(disty%in%c(4,6)|is.null(inits$disp))
      disp <- as.numeric(inits$disp)
    else
      disp <- rep(-999999,length(alpha))
    if(return_list) res <- list(alpha=alpha,beta=beta,eta=eta,disp=disp)
    else res <- c(alpha,beta,eta,disp)
  }

  if(is.numeric(inits)){
    start <- 0
    alpha <- inits[start + 1:S]
    start <- start + S
    beta <- inits[start + 1:((G*np))]
    start <- start + (G*np)
    eta <- inits[start + 1:(G - 1)]
    start <- start + (G-1)
    if(disty%in%c(4,6))
      disp <- inits[start + 1:S]
    else
      disp <- rep(-999999,S)

    if(return_list) res <- list(alpha=alpha,beta=beta,eta=eta,disp=disp)
    else res <- c(alpha,beta,eta,disp)

  }

  return(res)

}

"set_control_sam" <- function(control){
  if(is.null(control))control <- species_mix.control()
  control
}


"standardise.X" <- function (mat){
  X = scale(as.matrix(mat))
  dat.means = apply(as.matrix(mat), 2, mean, na.rm = TRUE)
  dat.sds = apply(as.matrix(mat), 2, sd, na.rm = TRUE)
  return(list(X = X, dat.means = dat.means, dat.sds = dat.sds))
}

"calc_info_crit_sam" <-  function(tmp) {
    k <- length(unlist(tmp$coefs))
    tmp$BIC <- -2 * tmp$logl + log(tmp$n) * k
    tmp$AIC <- -2 * tmp$logl + 2 * k
    return( tmp)
}

"calc_post_probs_sam" <- function( pis, logCondDens)  {
    logPostProbs <- log( pis) + logCondDens #would be better to be working with log(pis) previously but c'est la vie
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

# "check_distribution_clean_data_sam" <- function(arch_form, sp_form, dat, distribution){
#
#   # returns numeric 0, 1 or 2. if zero no species intercepts, if 1 speices intercepts, if 2 species model.
#   species_int_coefs <- check_species_formula(sp_form)
#   print(species_int_coefs)
#   if(species_int_coefs!=2){
#     mod_dat <- clean_data_sam(dat, arch_form, NULL, distribution)
#   } else {
#     mod_dat <- clean_data_sam(dat, arch_form, sp_form, distribution)
#   }
#
#   if(all(distribution=="bernoulli", species_int_coefs==0))fit_disty <- "bernoulli"
#   stop('Bernoulli distribution now requires species specific intercepts - look at "SpeciesMix" package for joint intercept estimation approach')
#   if(all(distribution=="bernoulli", species_int_coefs==1))fit_disty <- "bernoulli"
#   if(all(distribution=="bernoulli", species_int_coefs==2)){
#     fit_disty <- "bernoulli_partial"
#     stop('partial SAMs for a Bernoulli distribution has not been implemented yet - watch this space')
#   }
#
#   if(all(distribution=="poisson", species_int_coefs==0))
#     stop('Poisson distribution requires independent species intercepts')
#   if(all(distribution=="poisson", species_int_coefs==1))fit_disty <- "poisson"
#   if(all(distribution=="poisson", species_int_coefs==2)){
#     fit_disty <- "poisson_partial"
#     stop('partial SAMs for a Poisson distribution has not been implemented yet - watch this space')
#   }
#
#   if(all(distribution=="ippm", species_int_coefs==0))
#     stop('IPPM distribution requires independent species intercepts')
#   if(all(distribution=="ippm", species_int_coefs==1)) fit_disty <- "ippm"
#   if(all(distribution=="ippm", species_int_coefs==2)){
#     fit_disty <- "ippm_partial"
#     stop('partial SAMs for a IPPM distribution has not been implemented yet - watch this space')
#   }
#
#   if(all(distribution=="negative_binomial", species_int_coefs==0))
#     stop('Negative binomial distribution requires independent species intercepts')
#   if(all(distribution=="negative_binomial", species_int_coefs==1)) fit_disty <- "negative_binomial"
#   if(all(distribution=="negative_binomial", species_int_coefs==2)){
#     fit_disty <- "negative_binomial_partial"
#     stop('partial SAMs for a Negative Binomial distribution has not been implemented yet - watch this space')
# }
#
#   if(all(distribution=="tweedie", species_int_coefs==0))
#     stop('Tweedie distribution requires independent species intercepts')
#   if(all(distribution=="tweedie", species_int_coefs==1)) fit_disty <- "tweedie"
#   if(all(distribution=="tweedie", species_int_coefs==2)){
#     fit_disty <- "tweedie_partial"
#     stop('partial SAMs for a Tweedie distribution has not been implemented yet - watch this space')
#   }
#
#   return(list(fit_disty=fit_disty,mod_dat=mod_dat))
#
# }

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

"clean_data_sam" <- function(data, form1, form2, distribution){
    if(distribution=='ippm') na_rule <- "na.pass"
    else na_rule <- "na.exclude"
    mf.X <- stats::model.frame(form1, data = data, na.action = na_rule)
    if( !is.null( form2)){
      mf.W <- stats::model.frame(form2, data = data, na.action = na_rule)
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





