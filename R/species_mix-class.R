##### Species mix functions ####

#' @title species_mix objects
#' @rdname species_mix-class
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
#' @param archetype_fromula an object of class "formula" (or an object that can be
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
#' @param model_data a matrix of dataframe which contains the 'species_data'
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
#'  in the log-likelihood calculations. If NULL (default) then all weights are
#'  assumed to be identically 1. Because we are estimating the log-likelihood
#'  over species (rather than sites), the weights should be a vector n species
#'  long. The exception is under the use of the 'ippm' distribution where
#'  weights must be a nrow(data)*n_species matrix, which provides a
#'  species-specific background weights used to estimate the species-specific
#'  marginal likelihoods.
#' @param control a list of control parameters for optimisation and calculation.
#' See details. From \code{species_mix.control} for details on optimistaion
#' parameters.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficents. These are distribution specific, but at a
#' minimum you will need pis (additive_logitic transformed), alphas
#' (intercepts) and betas (mixing coefs).
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
#' library(ecomix)
#' set.seed(42)
#' sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),collapse = ','),")~1+x1+x2"))
#' sp_form <- ~ 1
#' theta <- matrix(c(1,-2.9,-3.6,1,-0.9,1,1,.9,7.9),3,3,byrow=TRUE)
#' dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),x2=stats::rnorm(100,0,2.5))
#' dat[,-1] <- scale(dat[,-1])
#' simulated_data <- simulate_species_mix_data(archetype_formula=sam_form, species_formula=sp_form,
#'                                             dat,theta,dist="bernoulli")
#' model_data <- make_mixture_data(species_data = simulated_data$species_data,
#'                                 covariate_data = simulated_data$covariate_data[,-1])
#' fm1 <- species_mix(sam_form, sp_form, model_data, distribution = 'bernoulli',
#'  n_mixtures=3)
#' }

"species_mix" <- function(archetype_formula = NULL, species_formula = stats::as.formula(~1), data,
                          n_mixtures = 3, distribution="bernoulli", offset=NULL,
                          weights=NULL, control=species_mix.control(), inits=NULL,
                          standardise = FALSE, titbits = TRUE){

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
  m <- match(c("data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  # if(distribution=="ippm")
  mf$na.action <- "na.pass"
  # else mf$na.action <- "na.exclude"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))

  # get the model matrix and find the fitting formula.
  # dist_dat <- check_distribution_clean_data_sam(archetype_formula, species_formula, mf, distribution)
  # fit_distribution <- dist_dat[[1]]
  # dat <- dist_dat[[2]]
  # print(fit_distribution)

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
  # disty.cases <- c("bernoulli","bernoulli_sp","poisson","ippm",
  #                  "negative_binomial","tweedie","gaussian")
  disty.cases <- c("bernoulli","poisson","ippm","negative_binomial","tweedie","gaussian")
  disty <- get_distribution_sam(disty.cases, distribution)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  weights <- get_weights_sam(mf,species_names,distribution)
  # cat(colnames(weights),"\n")

  if(distribution=='ippm'){
    if(!all(colnames(y)%in%colnames(weights)))
      # cat(colnames(y),"\n")
      # cat(colnames(weights),"\n")
      stop('When modelling a inhomogenous poisson point process model,
           weights colnames must match species data colnames')
    if(any(dim(y)!=dim(weights)))
      stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
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
  tmp <- species_mix.fit(y=y, X=X, W=W, G=n_mixtures, weights=weights, offset=offset,
                         distribution_numeric=disty, control, y_is_na=y_is_na, inits=inits,
                         archetype_formula=archetype_formula,species_formula=species_formula)

  tmp$dist <- disty.cases[disty]

  tmp$pis <- ecomix:::additive_logistic(tmp$eta)

  if(n_mixtures>1)
    tmp$post_probs <- ecomix:::calc_post_probs_sam(tmp$pis,tmp$loglikeSG)

  #Information criteria
  tmp <- calc_info_crit_sam(tmp)

  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_sam(titbits, y, X, W=NULL, weights, offset, archetype_formula, species_formula, control, disty.cases[disty])
  class(tmp) <- c("species_mix",distribution)
  return(tmp)
}

#'@rdname species_mix-class
#'@name species_mix.fit
#'@param y is a matrix genertated from \link[stats]{model.response} containing the species information. The matrix has the dimensions n_sites * n_species.
#'@param X is a design matrix for the archetype_formula dimension n_sites * n_covariates.
#'@param W is a design matrix for species_formula and will be implemented if species_formula has covariates.
#'@param G is the number of species archetypes that are being estimated.
#'@param weights is used in alternative way depending on the error distribution used. See \link[ecomix]{species_mix} for more details.
#'@param distribution_numeric the error distribution to used in species_mix estimation. Currently, 'bernoulli', 'poisson', 'ippm' (Poisson point process), 'negative_binomial' and 'tweedie' are avaliable - internal conversion of distribution to a integer.
#'@param offset this is a vector of site specific offsets, this might be something like area sampled at sites.
#'@param control this is a list of control parameters that alter the specifics of model fitting. See \link[ecomix]{species_mix.control} for details.
#'@param y_is_na This is a logical matrix used specifically with 'ippm' modelling - don't worry about this, it'll be worked out for you. Yay!
"species_mix.fit" <- function(y, X, W=NULL, G, weights, offset, distribution_numeric, control, y_is_na=NULL, inits=NULL,
                              archetype_formula=NULL,species_formula=NULL){

  disty.cases <- c("bernoulli","poisson","ippm",
                   "negative_binomial","tweedie","gaussian")
  distribution <- disty.cases[distribution_numeric]
  if(!any(distribution==c("bernoulli","poisson","ippm",
                          "negative_binomial","tweedie","gaussian")))
    stop('current only please check the distribution you are fitting')
  sp.form <- update(archetype_formula,obs~1+.)

  # need to insert two new calls:
  # get_starting_values_sam() which will generate starting values based on EM or clustering of coefs.
  # sam_optimise() which will fit the model in cpp
  S <- ncol(y) # how many species are there?

  starting_values <- get_start_vals_sam(y=y, X=X, weights=weights, offset=offset,
                                        disty=distribution_numeric,
                                        G=G, S=S, y_is_na=y_is_na, control=control)

  # spp_wts <- starting_values$spp_wts
  # site_spp_wts <- starting_values$site_spp_wts
  if(is.null(y_is_na))y_is_na <- matrix(FALSE,nrow(y),ncol(y)) #generate y_is_na for for c++ code.

  tmp <- sam_optimise(y, X, offset, starting_values$spp_wts, starting_values$site_spp_wts,
                      y_is_na, starting_values$nS, starting_values$nG, starting_values$nObs,
                      distribution_numeric, starting_values, control)

  # tmp <- switch(distribution,
  #               # bernoulli = species_mix_bernoulli(sp.form, y, X, G, inits,
  #                                                 # control),
  #               bernoulli = species_mix_bernoulli_sp(y = y, X = X, offset = offset,
  #                                                  weights = weights,
  #                                                  G = G, control = control),
  #               poisson = species_mix_ippm(y = y, X = X, weights = matrix(1,nrow(y),ncol(y)),
  #                                    offset = offset,  G = G, control = control,
  #                                    y_is_na = matrix(1,nrow(y),ncol(y))),
  #               ippm = species_mix_ippm(y = y, X = X, weights = weights,  offset = offset,
  #                                 G = G, control = control, y_is_na = y_is_na),
  #               negative_binomail = species_mix_nbinom(sp.form, y, X, G, inits,
  #                                                     control),
  #               tweedie = species_mix_tweedie(sp.form, y, X, G, inits,
  #                                            control),
  #               gaussian = species_mix_gaussian(sp.form, y, X, G, inits,
  #                                              controls))
  return(tmp)
}

#'@rdname species_mix-class
#'@name species_mix.fit
#'@examples
#' \dontrun{
#' fmods <- species_mix.multifit(sam_form, sp_form, model_data, distribution = 'bernoulli', nstarts = 10, n_mixtures=3)
#' }
"species_mix.multifit" <- function(archetype_formula = NULL, species_formula = ~1, data, distribution="bernoulli",
                                   nstart = 10, n_mixtures = 3, offset=NULL, weights=NULL,
                                   control=species_mix.control(), inits=NULL, standardise = FALSE){
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


  # Create model matrix
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  if(distribution=="ippm") mf$na.action <- "na.pass"
  else mf$na.action <- "na.exclude"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))

  # get the model matrix and find the fitting formula.
  dist_dat <- check_distribution_clean_data_sam(archetype_formula, species_formula, mf, distribution)
  fit_distribution <- dist_dat[[1]]
  dat <- dist_dat[[2]]

  # get responses
  y <- stats::model.response(dat$mf.X)

  # logical matirx needed for removing NAs from response and weights.
  if(distribution=='ippm')y_is_na <- is.na(y)
  else y_is_na <- NULL

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
  W <- get_W_sam(species_formula, dat$mf.W)

  #get distribution
  disty.cases <- c("bernoulli","bernoulli_sp","poisson","ippm",
                   "negative_binomial","tweedie","gaussian")
  print(fit_distribution)
  disty <- get_distribution_sam(disty.cases, fit_distribution)

  # get offsets and weights
  offset <- get_offset_sam(mf)
  weights <- get_weights_sam(mf,S,distribution)

  if(distribution=='ippm'){
    if(!all(colnames(y)%in%colnames(weights)))
      stop('When modelling a inhomogenous poisson point process model,
           weights colnames must match species data colnames')
    if(any(dim(y)!=dim(weights)))
      stop('When modelling a inhomogenous poisson point process model,
           weights needs to have the same dimensions at the
           species data - n_sites x n_species')
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

  tmp_fun <- function(x){
      if( !control$quiet & nstart>1)
        setTxtProgressBar(pb, x)
      tmpQuiet <- control$quiet
      control$quiet <- TRUE
      dumbOut <- capture.output(tmp <- species_mix.fit(y=y, X=X, W=W, weights=weights, offset=offset,
                                                       distribution_numeric=disty, n_mixtures=n_mixtures,
                                                       inits = inits, control=control, y_is_na=y_is_na,
                                                       estimate_variance=control$control$est_var))
      control$quiet <- tmpQuiet
      tmp$dist <- disty.cases[disty]

      #Information criteria
      tmp <- calc_info_crit_sam(tmp)

      #titbits object, if wanted/needed.
      tmp$titbits <- get_titbits_sam(titbits, y, X, W, weights, offset, archetype_formula, species_formula, control, disty.cases[disty])
      # tmp$titbits$disty <- disty

      #the last bit of the regional_mix object puzzle
      # tmp$call <- call
      class(tmp) <- c("species_mix", distribution)
      return( tmp)
   }

   #Fit the model many times
   many_starts <- surveillance::plapply(seq_len(nstart), tmp_fun, .parallel = control$cores, .verbose = !control$quiet)
   return(many_starts)
}

"get_start_vals_sam" <- function(y, X, weights, offset, disty, G, S, y_is_na, control){

  if(!control$quiet)
    print( "Obtaining starting values...")

  if(disty==1){
    tmp <- get_starting_values_bernoulli_sp(y=y, X=X, weights=weights, offset=offset,
                                            G=G, S=S, control=control)
  }
  if(disty==2){
    tmp <- get_starting_values_poisson(y=y, X=X, weights=weights, offset=offset,
                                       G=G, S=S, control=control)
  }
  if(disty==3){
    tmp <- get_starting_values_ippm(y=y, X=X, weights=weights, offset=offset,
                                    G=G, S=S, y_is_na=y_is_na, control=control)
  }
  if(disty==4){
    tmp <- get_intital_values_negative_binomial()
  }
  if(disty==5){
    tmp <- get_initial_values_tweedie()
  }
  if(disty==6){
    tmp <- get_initial_values_gaussian()
  }

  return(tmp)
}


"sam_optimise" <- function(y, X, offset, spp_wts, site_spp_wts, y_is_na, nS, nG, nObs, disty, start_vals, control){

  inits <- c(start_vals$alphas, start_vals$betas, start_vals$eta, start_vals$disp)
  np <- as.integer(ncol(X[,-1]))
  n <- as.integer(nrow(X))

  # parameters to optimise
  alpha <- as.numeric(start_vals$alphas)
  beta <- as.numeric(start_vals$betas)
  eta <- as.numeric(start_vals$eta)
  disp <- as.numeric(start_vals$disp)

  #scores
  alpha.score <- as.numeric(rep(NA, length(alpha)))
  beta.score <- as.numeric(rep(NA, length(beta)))
  eta.score <- as.numeric(rep(NA, length(eta)))
  disp.score <- as.numeric(rep(NA, length(disp)))
  getscores <- 1
  scores <- as.numeric(rep(NA,length(c(alpha,beta,eta,disp))))

  #model quantities
  pis_out <- as.numeric(rep(NA, nG))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c( nObs, nS, nG)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, nS))
  loglikeSG  <- as.numeric(matrix(NA, nrow = nS, ncol = nG))

  #c++ call to optimise the model (needs pretty good starting values)
  tmp <- .Call("species_mix_cpp",
               as.numeric(as.matrix(y)), as.numeric(as.matrix(X[,-1])), as.numeric(offset), as.numeric(spp_wts),
               as.numeric(as.matrix(site_spp_wts)), as.integer(as.matrix(!y_is_na)),
               # SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_wts, SEXP Rsite_spp_wts, SEXP Ry_not_na, // data
               as.integer(nS), as.integer(nG), as.integer(np), as.integer(nObs), as.integer(disty),
               # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
               as.double(alpha), as.double(beta), as.double(eta), as.double(disp),
               # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
               alpha.score, beta.score, eta.score, disp.score, as.integer(control$getscores), scores,
               # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
               pis_out, mus, loglikeS, loglikeSG,
               # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
               as.integer(control$maxit_cpp), as.integer(control$trace_cpp), as.integer(control$nreport_cpp),
               as.numeric(control$abstol_cpp), as.numeric(control$reltol_cpp), as.integer(control$conv_cpp), as.integer(control$printparams_cpp),
               # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
               as.integer( control$optimise_cpp), as.integer(control$loglOnly_cpp), as.integer( control$derivOnly_cpp), as.integer(control$optiDisp),
               # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
               PACKAGE = "ecomix")

  ret <- tmp
  ret$logl <- ret$logl * -1
  ret$mus <- array(mus, dim=c(nObs, nS, nG))
  ret$coefs <- list(alpha = ret$alpha, beta = ret$beta, eta = ret$eta, disp = ret$disp)
  ret$scores <- list(alpha.scores = alpha.score, beta.scores = beta.score, eta.scores=eta.score, disp.scores=disp.score)
  ret$S <- nS; ret$G <- nG; ret$np <- np; ret$n <- nObs;
  ret$start.vals <- inits
  ret$loglikeSG <- matrix(loglikeSG,  nrow = nS, ncol = nG)  #for residuals
  ret$loglikeS <- loglikeS  #for residuals
  return(ret)
}

#'@rdname species_mix-class
#'@name control
#'@param quite Should any reporting be performed? Default is FALSE, for reporting.
#'@param trace int 1=model will report parameter estimates and loglikelihood at each iteration. 0=quite.
#'@param reltol function that determines the relative tolernace for model convergence. Default is quite strict.
#'@param maxit Maximum number of evaluations of the objective function allowed. Defaults to 500.
#'@param cores The number of cores to use in fitting of species mix models. These will be largely used to model the species-specific parameteres.
#'@param em_prefit Logical if TRUE the model will run a slower EM algorithim fit to find starting values.
#'@param em_steps int Default is 3, the number of EM iterations to get to starting values.
#'@param em_refit int Default is 1, number of times to refit using EM.
#'@param calculate_hessian logical if TRUE model will numerically estimate the variance covariance matrix.
#'@param residuals logical if TRUE model will estimate residuals.
#'@export
"species_mix.control" <- function(maxit = 1000,
                                  quiet = FALSE,
                                  trace = 1,
                                  cores = 1,
                                  residuals = FALSE,
                                  ## intialisation controls
                                  init_method = 'kmeans',
                                  init_sd = 1,
                                  ## EM algorithim controls
                                  em_prefit = TRUE,
                                  em_steps = 3,
                                  em_refit = 3,
                                  em_maxit = 3,
                                  em_abstol = sqrt(.Machine$double.eps),
                                  em_reltol = reltol_fun,
                                  em_maxtau = 0.8,
                                  em_calculate_hessian = FALSE,
                                  em_full_model = FALSE,
                                  r1=1,
                                  ## c++ controls
                                  maxit_cpp = 1000,
                                  trace_cpp = 0,
                                  nreport_cpp = 1,
                                  abstol_cpp = sqrt(.Machine$double.eps),
                                  reltol_cpp = sqrt(.Machine$double.eps),
                                  conv_cpp = 1,
                                  printparams_cpp = 0,
                                  optimise_cpp = 1,
                                  loglOnly_cpp = 0,
                                  derivOnly_cpp = 0,
                                  getscores_cpp = 0,
                                  calculate_hessian_cpp = TRUE,
  ...){
               #general controls
  rval <- list(maxit = maxit, quiet = quiet, trace = trace,
               cores = cores,  residuals = residuals,
               #initialisation controls
               init_method = init_method, init_sd = init_sd,
               #em controls
               em_prefit = em_prefit, em_refit = em_refit, em_steps = em_steps, em_maxit = em_maxit,
               em_abstol = em_abstol,
               em_reltol = em_reltol, em_maxtau = em_maxtau, em_calculate_hessian = em_calculate_hessian,
               em_full_model = em_full_model, r1 = r1,
               #cpp controls
               maxit_cpp = maxit_cpp, trace_cpp = trace_cpp, nreport_cpp = nreport_cpp,
               abstol_cpp = abstol_cpp, reltol_cpp = reltol_cpp, conv_cpp = conv_cpp,
               printparams_cpp = printparams_cpp, optimise_cpp = optimise_cpp,
               loglOnly_cpp = loglOnly_cpp, derivOnly_cpp = derivOnly_cpp,
               getscores_cpp = getscores_cpp, calculate_hessian_cpp = calculate_hessian_cpp)
  rval <- c(rval, list(...))
  if (is.null(rval$em_reltol))
    rval$em_reltol <- sqrt(.Machine$double.eps)
  rval
}

#' @rdname species_mix-class
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

    sp_weights <- lapply(seq_along(sp_name),function(x)(weights=df$area/as.numeric(species_specific_cell_counts[[x]][match(df$id,as.numeric(names(species_specific_cell_counts[[x]])))])))

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
    if (distribution == "tweedie") {
      tmp <- rep(6e+05, dim(X)[1])
      while (max(tmp, na.rm = TRUE) > 5e+05 | sum(tmp) < 100) {
        theta[g, 1] <- stats::runif(1, -15, 5)
        theta[g, 1] <- stats::runif(1, 1, 5)
        sp.int[s] <- theta[g, 1]
        lgtp <- X %*% theta[g, ]
        p <- exp(lgtp)
        tmp <- rTweedie(dim(X)[1], mu = p, phi = 2, p = 1.6)
      }
      out[, s] <- tmp
    }
    if (distribution == "gaussian") {
      sp.int <- NULL
      tmp <- rep(1e+05, dim(X)[1])
      while (max(tmp, na.rm = TRUE) > 50000 | sum(tmp) < 100) {
        theta[g, 1] <- stats::runif(1, 1, 500)
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

#'@rdname species_mix-class
#'@name species_mix_estimate_groups
#'@description This function runs the 'species_mix' function but it iterates through groups (1 to 10) by default.
"species_mix_estimate_groups" <- function(archetype_formula = NULL, species_formula = ~1,
                                          data, n_mixtures = 1:10, distribution="bernoulli",
                                          offset=NULL, weights=NULL, control=species_mix.control(),
                                          inits=NULL, standardise = FALSE){
  my.fun <- function(G,archetype_formula,species_formula,data,offset,weights,distribution,inits,control){
    message("Fitting group",G,"\n")

    ecomix::species_mix(archetype_formula, species_formula, data, n_mixtures)

  }
  out <- surveillance::plapply(G, my.fun, form, dat, .parallel = control$cores, .verbose = !control$quiet)
  aic <- rep(0,length(G))
  bic <- rep(0,length(G))
  fm <- list()
  for(i in seq_along(G))
    if(!is.atomic(out[[i]])){
      aic[i] <- out[[i]]$aic
      bic[i] <- out[[i]]$bic
      fm[[i]] <- list(logl=out[[i]]$logl,coef=out[[i]]$coef,tau=out[[i]]$tau,pi=out[[i]]$pi,covar=out[[i]]$covar)
    }
  return(list(aic=aic,bic=bic,fm=fm))
}

#'@rdname species_mix-class
#'@name species_mix.predict
#'@param object is a matrix model returned from the species_mix model.
#'@param new_obs a matrix of new observations for prediction.
#'@description Predicts SAM probabilities at a series of sites. Confidence intervals can be calculated if variance-covariance matrix is estimated during species_mix model fit.
#'@examples
#'\dontrun{
#'fm1 <- species_mix(form,data)
#'preds_fm1 <- predict(fm1,newdata)}
"species_mix.predict" <-function (object, new_obs, ...){
  mixture.model <- object
  if (class(mixture.model)[2] == "bernoulli") {
    G <- length(mixture.model$pi)
    covar <- mixture.model$covar[-(1:(G - 1)), -(1:(G -
        1))]
    coef <- mixture.model$coef
    model.fm <- stats::as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- stats::model.matrix(model.fm, new_obs)
    link.fun <- stats::make.link("logit")
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
    model.fm <- stats::as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- cbind(stats::model.matrix(model.fm, new_obs), 1)
    offset <- stats::model.frame(model.fm, data = new_obs)
    offset <- stats::model.offset(offset)
    if (is.null(offset))
      offset <- rep(0, nrow(X))
    outvar <- matrix(NA, dim(X)[1], G)
    outpred <- matrix(NA, dim(X)[1], G)
    colnames(outvar) <- colnames(outpred) <- paste("G",
      1:G, sep = ".")
    for (g in 1:G) {
      s.outvar <- matrix(NA, dim(X)[1], length(sp.int))
      s.outpred <- matrix(NA, dim(X)[1], length(sp.int))
      for (s in seq_along(sp.int)) {
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
  if (class(mixture.model)[2] == "ippm" | class(mixture.model)[2] == "poisson") {
    G <- length(mixture.model$pi)
    covar <- mixture.model$covar[-1 * c(1:(G - 1)), -1 * c(1:(G - 1))]
    sp.int <- mixture.model$sp_intercept
    coef <- mixture.model$coef
    model.fm <- stats::as.formula(mixture.model$formula)
    model.fm[[2]] <- NULL
    X <- cbind(stats::model.matrix(model.fm, new_obs))
    offset <- stats::model.frame(model.fm, data = new_obs)
    offset <- stats::model.offset(offset)
    if (is.null(offset))
      offset <- rep(0, nrow(X))
    outvar <- matrix(NA, dim(X)[1], G)
    outpred <- matrix(NA, dim(X)[1], G)
    colnames(outvar) <- colnames(outpred) <- paste("G", 1:G, sep = ".")
    for (g in 1:G) {
      s.outvar <- matrix(NA, dim(X)[1], length(sp.int))
      s.outpred <- matrix(NA, dim(X)[1], length(sp.int))
      for (s in seq_along(sp.int)) {
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

###### SAM internal functions #####

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
# "lambda_penalisation_fun" <- function(x,kappa=0.1){
#   res <- min(x,na.rm = TRUE)+kappa*(max(x,na.rm = TRUE)-min(x,na.rm = TRUE))
#   res
# }

"lambda_penalisation_fun" <- function(x,lambda,kappa=0.1){ #assumes that x spans to pretty-well the unpenalised estiamtes
  min.effective.penalty <- min( which( abs( x-tail( x, 1)) < 0.01 * abs( tail( x, 1))))    #the first that lambda that gives a coef close to the last lambda's corresponding coef
  min.effective.penalty <- lambda[min.effective.penalty]
  target.penalty <- kappa * min.effective.penalty
  res.pos <- which.min( (lambda-target.penalty)^2)
  res <- x[res.pos]
  return( res)
}

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

"get_titbits_sam" <- function( titbits, y, X, W, weights, offset,
                               archetype_formula, species_formula, control,
                               distribution)  {
    if( titbits==TRUE)
      titbits <- list( Y = y, X = X, W = W, weights=weights, offset=offset,
                       archetype_formula =  archetype_formula,
                       species_formula = species_formula,
                       control = control,
                       distribution = distribution)
    else{
      titbits <- list()
      if( "Y" %in% titbits)
        titbits$Y <- y
      if( "X" %in% titbits)
        titbits$X <- X
      if( "W" %in% titbits)
        titbits$W <- W
      if( "offset" %in% titbits)
        titbits$offset <- offset
      if( "weights" %in% titbits)
        titbits$weights <- weights
      if( "archetype_formula" %in% titbits)
        titbits$archetype_formula <- archetype_formula
      if( "species_formula" %in% titbits)
        titbits$species_formula <- species_formula
      if( "control" %in% titbits)
        titbits$control <- control
      if( "distribution" %in% titbits)
        titbits$distribution <- distribution
    }
    return( titbits)
}

"get_weights_sam"  <- function(mf,sp_names,distribution){
  if(distribution=='ippm'){
    weights <- model.weights(mf)
    if(!is.null(weights)){
      weights <- subset(weights, select = colnames(weights)%in%sp_names)
    } else {
      weights <- matrix(1,nrow(mf),length(sp_names))
    }
  } else {
    weights <- rep(1, nrow(mf))
  }
  return(weights)
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

"reltol_fun" <- function(logl_n1, logl_n){
  return(abs(logl_n1 - logl_n) > (abs(logl_n1 - logl_n) / abs(logl_n)))
}

"set_control_sam" <- function(control){
  if(is.null(control))control <- species_mix.control()
  control
}

# "clean_data_sam" <- function(form, data){
#
#   mf.X <- stats::model.frame(formula = form, data = data, na.action = na.exclude)
#   ids <- rownames(mf.X)
#   res <- list(ids=ids, mf.X=mf.X)
#   return(res)
# }

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

"check_distribution_clean_data_sam" <- function(arch_form, sp_form, dat, distribution){

  # returns numeric 0, 1 or 2. if zero no species intercepts, if 1 speices intercepts, if 2 species model.
  species_int_coefs <- check_species_formula(sp_form)
  print(species_int_coefs)
  if(species_int_coefs!=2){
    mod_dat <- clean_data_sam(dat, arch_form, NULL, distribution)
  } else {
    mod_dat <- clean_data_sam(dat, arch_form, sp_form, distribution)
  }

  if(all(distribution=="bernoulli", species_int_coefs==0))fit_disty <- "bernoulli"
  stop('Bernoulli distribution now requires species specific intercepts - look at "SpeciesMix" package for joint intercept estimation approach')
  if(all(distribution=="bernoulli", species_int_coefs==1))fit_disty <- "bernoulli"
  if(all(distribution=="bernoulli", species_int_coefs==2)){
    fit_disty <- "bernoulli_partial"
    stop('partial SAMs for a Bernoulli distribution has not been implemented yet - watch this space')
  }

  if(all(distribution=="poisson", species_int_coefs==0))
    stop('Poisson distribution requires independent species intercepts')
  if(all(distribution=="poisson", species_int_coefs==1))fit_disty <- "poisson"
  if(all(distribution=="poisson", species_int_coefs==2)){
    fit_disty <- "poisson_partial"
    stop('partial SAMs for a Poisson distribution has not been implemented yet - watch this space')
  }

  if(all(distribution=="ippm", species_int_coefs==0))
    stop('IPPM distribution requires independent species intercepts')
  if(all(distribution=="ippm", species_int_coefs==1)) fit_disty <- "ippm"
  if(all(distribution=="ippm", species_int_coefs==2)){
    fit_disty <- "ippm_partial"
    stop('partial SAMs for a IPPM distribution has not been implemented yet - watch this space')
  }

  if(all(distribution=="negative_binomial", species_int_coefs==0))
    stop('Negative binomial distribution requires independent species intercepts')
  if(all(distribution=="negative_binomial", species_int_coefs==1)) fit_disty <- "negative_binomial"
  if(all(distribution=="negative_binomial", species_int_coefs==2)){
    fit_disty <- "negative_binomial_partial"
    stop('partial SAMs for a Negative Binomial distribution has not been implemented yet - watch this space')
}

  if(all(distribution=="tweedie", species_int_coefs==0))
    stop('Tweedie distribution requires independent species intercepts')
  if(all(distribution=="tweedie", species_int_coefs==1)) fit_disty <- "tweedie"
  if(all(distribution=="tweedie", species_int_coefs==2)){
    fit_disty <- "tweedie_partial"
    stop('partial SAMs for a Tweedie distribution has not been implemented yet - watch this space')
  }

  return(list(fit_disty=fit_disty,mod_dat=mod_dat))

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

# "apply_glm" <-  function (i,form,datsp,tau,n){
#     dat.tau <- rep(tau[,i],each=n)
#     x <- stats::model.matrix(stats::as.formula(form),data=datsp)
#     y <- datsp$obs
#     f.mix <- stats::glm.fit(x,y,dat.tau,family=stats::binomial())
#     return(list(coef=f.mix$coef))
#   }


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


# "create_starting_values" <- function (S,G,n,form,datsp,control){
# 	environment(form) <- environment()
#     tau <- matrix(stats::runif(S*G),S,G)
#     tau <- (tau/rowSums(tau))
#     fmM <- list()
#     for(i in 1:G){
#       pi[i] <- sum(tau[,i])/S
#     }
#     fmM <- surveillance::plapply(1:G,apply_glm,form,datsp,tau,n,.parallel=control$cores, .verbose = !control$quiet)
#     first.fit <- list(x=stats::model.matrix(stats::as.formula(form),data=datsp),y=datsp$obs,formula=form)
#     return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
#   }



# "distr.binom" <- function( p){
#     nobs <- length( p)
#     new.dist <- old.dist <- rep( 0, nobs+1)
#     old.dist[1] <- 1-p[1]
#     old.dist[2] <- p[1]
#     for( ii in 2:nobs){
#       new.dist[1] <- old.dist[1]*(1-p[ii])
#       for( jj in 2:ii)
#         new.dist[jj] <- old.dist[jj-1]*p[ii] + old.dist[jj]*(1-p[ii])
#       new.dist[ii+1] <- old.dist[ii]*p[ii]
#       old.dist <- new.dist
#     }
#     return( new.dist)
#   }


# "dPoisGam" <- function ( y, lambda, mu.Z, alpha, LOG=TRUE){
#     #function to calculate Random sum (Tweedie) densities.
#     #y is the value of the r.v.  Can be a vector
#     #mu.N is the mean of the Poisson summing r.v. Can be a vector of length(y)
#     #mu.Z is the mean of the Gamma rv Can be a vector of length(y)
#     #alpha is the `other' parameter of the gamma distribution s.t. var = ( mu.Z^2)/alpha Can be a vector of length(y)
#     #If mu.N, mu.Z or alpha are scalare but y isn't then they will be used for all y. If lengths mis-match then error
#     #LOG=TRUE gives the density on the log scale
#     #do.checks=TRUE checks the input vectors for compatability and gives errors / changes them as appropriate.
#     #do.checks=FALSE doesn't check and relies on the user to have things right. If not right then catastrophic failure may occur.
#     mu.N <- lambda
#     if( !all( is.element( c( length( mu.N), length( mu.Z), length( alpha)), c( length( y), 1)))){
#       print( "Error: length of parameter vectors does not match length of random variable vector")
#       print( "Exitting")
#       return()
#     }
#
#     if( length( mu.N) != length( y))
#       mu.N <- rep( mu.N, length( y))
#     if( length( mu.Z) != length( y))
#       mu.Z <- rep( mu.Z, length( y))
#     if( length( alpha) != length( y))
#       alpha <- rep( alpha, length( y))
#
#     res <- .Call( "dTweedie", as.numeric( y), as.numeric( mu.N), as.numeric( mu.Z), as.numeric( alpha), as.integer( LOG),PACKAGE="ecomix")
#
#     return( res)
#
#   }


# "dPoisGamDerivs" <- function ( y=NULL, lambda=NULL, mu.Z=NULL, alpha=NULL, do.checks=TRUE){
#     #function to calculate Random sum (Tweedie) densities.
#     #y is the value of the r.v.  Can be a vector
#     #mu.N is the mean of the Poisson summing r.v. Can be a vector of length(y)
#     #mu.Z is the mean of the Gamma rv Can be a vector of length(y)
#     #alpha is the `other' parameter of the gamma distribution s.t. var = ( mu.Z^2)/alpha Can be a vector of length(y)
#     #If mu.N, mu.Z or alpha are scalare but y isn't then they will be used for all y. If lengths mis-match then error
#     #LOG=TRUE gives the density on the log scale
#     #do.checks=TRUE checks the input vectors for compatability and gives errors / changes them as appropriate.
#     #do.checks=FALSE doesn't check and relies on the user to have things right. If not right then catastrophic failure may occur.
#
#     mu.N <- lambda
#     if( do.checks){
#       if( any( is.null( c( y, mu.N, mu.Z, alpha)))){
#         print( "Error: null input values -- please check.  Null values are:")
#         tmp <- double( is.null( c( y, mu.N, mu.Z, alpha)))
#         names( tmp) <- c( "y", "mu.N","mu.Z","alpha")
#         print( tmp)
#         print( "Exitting")
#         return()
#       }
#
#       if( !all( is.element( c( length( mu.N), length( mu.Z), length( alpha)), c( length( y), 1)))){
#         print( "Error: length of parameter vectors does not match length of random variable vector")
#         print( "Exitting")
#       }
#
#       if( length( mu.N) != length( y))
#         mu.N <- rep( mu.N, length( y))
#       if( length( mu.Z) != length( y))
#         mu.Z <- rep( mu.Z, length( y))
#       if( length( alpha) != length( y))
#         alpha <- rep( alpha, length( y))
#     }
#
#     res <- .Call( "dTweedieDeriv", as.numeric( y), as.numeric( mu.N), as.numeric( mu.Z), as.numeric( alpha),PACKAGE="ecomix")
#     colnames( res) <- c("lambda","mu.Z","alpha")
#     return( res)
#
#   }

#
# "dTweedie" <- function ( y, mu, phi, p, LOG=TRUE){
#     lambda <- ( mu^( 2-p)) / ( phi*(2-p))
#     alpha <- ( 2-p) / ( p-1)
#     tau <- phi*(p-1)*mu^(p-1)
#     mu.Z <- alpha * tau
#
#     dens <- dPoisGam( y, lambda, mu.Z, alpha, LOG)
#
#     return( dens)
#
#   }


# "estimate_pi" <- function (j,sp,spname,datsp,fmM,pi,G,first.fit){
#
#     tmp.like <- rep(0,G)
#     tau <- rep(0,G)
#     sel.sp <- which(sp==spname[j])
#
#     link.fun <- stats::make.link("logit")
#     for(i in 1:G) {
#       if(length(fmM[[i]]$coef)==1){
#         lpre <- link.fun$linkinv(first.fit$x[sel.sp,]*fmM[[i]]$coef)
#       }else{    lpre <- link.fun$linkinv(first.fit$x[sel.sp,]%*%fmM[[i]]$coef)}
#       obs <- datsp$obs[sel.sp]
#       lpre[obs==0]<- 1- lpre[obs==0]
#       tmp.like[i] <- sum(log(lpre))
#     }
#     eps <- max(tmp.like)
#     sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
#     for(i in 1:G) {
#       tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
#     }
#
#     return(list(tau=tau,sum.like=sum.like))
#   }



# "fitmix" <-  function (form,datsp,sp,G=2,control) {
#     ## dat2 has colums obs,sp
#     ##
#     temp.warn <- getOption( "warn")
#     options( warn=-1)
#
#     sp.name <- unique(sp)
#     S <- length(unique(sp))
#     n <- length(which(sp==sp.name[1]))
#
#     if(!control$quiet){
#       if(control$trace==1){
#         message("Fitting Group",G,"\n")
#         message("Iteration | LogL \n")
#       }
#     }
#
#     dat.tau <- 0
#     pi <- rep(0,G)
#     ite <- 1
#     logL <- -99999999
#     old.logL <- -88888888
#
#     ## set up initial GLM
#     t1 <- create_starting_values(S,G,n,form,datsp,control)
#     pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
#
#     while(abs(logL-old.logL) > control$em_abstol & ite<=control$maxit){
#       old.logL <- logL
#
#       for(i in 1:G){
#         pi[i] <- sum(tau[,i])/S
#       }
#
#       if(any(pi==0)) {
#         cat("pi has gone to zero - restarting fitting \n")
#         t1 <- create_starting_values(S,G,n,form,datsp,control)
#         pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
#         ite <- 1
#       }
#
#       fmM <- surveillance::plapply(1:G,weighted_glm,first.fit,tau,n,fmM,sp,.parallel=control$cores, .verbose = !control$quiet)
#
#
#       logL <- 0
#       tmp.like <- matrix(0,S,G)
#
#       est.tau <- surveillance::plapply(1:S,estimate_pi,sp,sp.name,datsp,fmM,pi,G,first.fit,.parallel=control$cores, .verbose = !control$quiet)
#
#       for(j in 1:S){
#         if(is.atomic(est.tau[[j]])){ print (est.tau[[j]])} else
#         {
#           tau[j,] <- est.tau[[j]]$tau
#           logL <- logL+est.tau[[j]]$sum.like
#         }
#       }
#
#       if(!control$quiet){
#       if(control$trace==1) cat(ite,"     | ",logL,"\n")
#       }
#       ite <- ite+1
#     }
#     fm.out <- data.frame(matrix(0,G,length(fmM[[1]]$coef)))
#     names(fm.out) <- names(fmM[[1]]$coef)
#     tau <- data.frame(tau)
#     names(tau) <- paste("grp.",1:G,sep="")
#     EN <- -sum(unlist(tau)*log(unlist(tau)))
#     d <- length(unlist(fm.out)) + length(tau)-1
#     for(i in 1:G) {
#       fm.out[i,] <- fmM[[i]]$coef
#       ##dat.tau[,i] <- rep(tau[,i],each=n)
#     }
#
#     names(pi) <- paste("G",1:G,sep=".")
#     t.pi <- additive_logistic(pi,TRUE)
#     parms <- c(t.pi[1:(G-1)],unlist(fm.out))
#     logL.full <- logL
#     logL <- logLmix(parms,first.fit,G,S,sp,sp.name)
#
#     options(warn=temp.warn)
#     if(control$em_full_model)  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,fmM=fmM,model.tau=dat.tau,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
#
#     return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))
#
#   }


# "fitmix.cpp" <- function (form,datsp,sp,G=2,pars=NA,control){
#     if(!is.numeric(sp)){
#       sp <- as.integer(factor(sp))
#     }
#     sp.name <- unique(sp)
#     S <- length(unique(sp))
#     n <- length(which(sp==sp.name[1]))
#     X <- stats::model.matrix(form, data = datsp[sp==sp.name[1],])
#     y <- datsp$obs
#     if(is.na(pars[1])) {
#       pars <- stats::runif(G-1+(ncol(X)*G),-2,2)
#       pars[1:(G-1)] <- stats::runif(G-1)
#     }
#     offset<-rep(0,n)
#     gradient <- rep(0,length(pars))
#     tau <- matrix(0,S,G) ##must leave this in as defines S & G
#     model.type<- as.integer(1)
#     loglike <- try(.Call("SpeciesMix",pars,y,X,sp,tau,gradient,offset,model.type,PACKAGE="ecomix"))
#
#     calc_deriv <- function(p){
#       gradient <- rep(0,length(pars))
#       ll <- .Call("Calculate_Gradient",p,y,X,sp,tau,gradient,offset,as.integer(1),PACKAGE="ecomix")
#       return(gradient)
#     }
#     hes <- 0
#     if(control$calculate_hessian_cpp){
#       hes <- numDeriv::jacobian(calc_deriv,pars)
#       dim(hes) <- rep(length(pars),2)
#       dim(hes) <- rep(length(pars),2)
#       rownames(hes) <- colnames(hes) <- c(paste("G.",1:(G-1),sep=""),paste("G",1:G,rep(colnames(X),each=G),sep="."))
#     }
#     if(!is.numeric(loglike)) loglike <- 0
#     pi <- pars[1:(G-1)]
#     coef <- pars[ (G):length(pars)]
#
#     r.logl <- logLmix(pars,list(y=y,x=stats::model.matrix(form, data = datsp)),G,S,sp,sp.name,out.tau=TRUE)
#     pi <- additive_logistic(pi)
#     names(pi) <- paste("G.",1:G,sep="")
#     coef <- matrix(coef,G,ncol(X))
#     rownames(coef) <- paste("G.",1:G,sep="")
#     colnames(coef) <- colnames(X)
#
#     list(logl=loglike,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4), n = length(pars), hessian=hes,gradient=gradient)
#   }




#
# "glm_fit_nbinom" <- function (x,y,offset=NULL,weights=NULL,mustart=NULL,control){
#     X <- x
#     if(is.null(offset)) offset <- 0
#     if(is.null(weights)) weights <- rep(1,length(y))
#     gradient <- rep(0,ncol(X)+1)
#     if(is.null(mustart)){ pars <- gradient+1} else{pars <- mustart}
#     fitted.values <- rep(0,length(y))
#     logl <- .Call("Neg_Bin",pars,X,y,weights,offset,gradient,fitted.values,PACKAGE="ecomix")
#     vcov <- 0
#     se <- rep(0,length(pars))
#     if(control$calculate_hessian_cpp) {
#       calc_deriv <- function(p){
#         gradient <- rep(0,length(pars))
#         ll <- .Call("Neg_Bin_Gradient",p,X,y,weights,offset,gradient,PACKAGE="ecomix")
#         return(gradient)
#       }
#       hes <- numDeriv::jacobian(calc_deriv,pars)
#       dim(hes) <- rep(length(pars),2)
#       vcov <- try(solve(hes))
#       se <- try(sqrt(diag(vcov)))
#       colnames(vcov) <- rownames(vcov) <- c("theta",colnames(X))
#     }
#     names(pars) <- names(se) <- names(gradient) <- c("theta",colnames(X))
#
#     return(list(logl=logl,coef=pars[-1],theta=pars[1],se=se[-1],se.theta=se[1],fitted=fitted.values,gradient=gradient,vcov=vcov))
#   }


# "nbglm" <-
#   function (form,data,weights=NULL,mustart=NULL,control)
#   {
#     X <- stats::model.matrix(form,data)
#     t1 <- stats::model.frame(form,data)
#     y <- stats::model.response(t1)
#     offset <- stats::model.offset(t1)
#     if(is.null(offset)) offset <- 0
#     if(is.null(weights)) weights <- rep(1,length(y))
#     gradient <- rep(0,ncol(X)+1)
#     if(is.null(mustart)){ pars <- gradient+1}else{pars <- mustart}
#     fitted.values <- rep(0,length(y))
#
#     logl <- .Call("Neg_Bin",pars,X,y,weights,offset,gradient,fitted.values,PACKAGE="ecomix")
#     vcov <- 0
#     se <- rep(0,length(pars))
#     if(control$calculate_hessian_cpp) {
#       calc_deriv <- function(p){
#         gradient <- rep(0,length(pars))
#         ll <- .Call("Neg_Bin_Gradient",p,X,y,weights,offset,gradient,PACKAGE="ecomix")
#         return(gradient)
#       }
#       hes <- numDeriv::jacobian(calc_deriv,pars)
#       dim(hes) <- rep(length(pars),2)
#       vcov <- try(solve(hes))
#       se <- try(sqrt(diag(vcov)))
#       colnames(vcov) <- rownames(vcov) <- c("theta",colnames(X))
#     }
#     names(pars) <- names(se) <- names(gradient) <- c("theta",colnames(X))
#
#     return(list(logl=logl,coef=pars[-1],theta=pars[1],se=se[-1],se.theta=se[1],fitted=fitted.values,gradient=gradient,vcov=vcov))
#   }


# "logLike.pars" <-
#   function (pi,coef,sp.form,sp.data,covar.data)
#   {
#
#     G <- length(pi)
#     pars <- c(additive_logistic(pi,TRUE)[1:(G-1)],unlist(coef))
#
#     S <- dim(sp.data)[2]
#     if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
#     n <- dim(sp.data)[1]
#     sp <- rep(sp.name,each=n)
#     var.names <- colnames(covar.data)
#     covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
#     names(covar.data) <- var.names
#     sp.data <- as.data.frame(sp.data)
#     data <- data.frame(obs=unlist(sp.data),covar.data)
#     names(data)[1] <- as.character(sp.form)[2]
#     first.fit <- list(y=data[,1],x=stats::model.matrix(sp.form,data=data))
#     logl <- -logLmix(pars,first.fit,G,S,sp,sp.name)
#     logl
#   }

# "logLmix" <-
#   function (pars,first.fit,G,S,sp,spname,out.tau=FALSE)
#   {
#     tau <- matrix(0,S,G)
#     ##tau,out.tau=FALSE
#     if(G>1){
#       fm <- pars[-1*(1:(G-1))]
#       pi <- pars[(1:(G-1))]
#       ##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
#       #pi <- tau[((length(tau)-(G-2)):length(tau))]
#       dim(fm) <- c(G,(length(pars)-(G-1))/G)
#       ##pi[G] <- 1-sum(pi)
#       pi <- additive_logistic(pi)
#     } else{
#       fm <- tau[1:(length(pars)-1)]
#       dim(fm) <- c(1,length(fm))
#       pi <- 1
#     }
#
#     link.fun <- stats::make.link("logit")
#
#
#     log.like <- 0
#     for(j in 1:S){
#       sel.sp <- which(sp==spname[j])
#       tmp.like <- rep(0,G)
#       for(i in 1:G){
#         if(length(fm[i,])==1){lpre <- first.fit$x[sel.sp,]*fm[i,]
#         }else{      lpre <- first.fit$x[sel.sp,]%*%fm[i,]}
#         tmp.like[i] <- sum(dbinom(first.fit$y[sel.sp],1,link.fun$linkinv(lpre),log=TRUE))
#       }
#       eps <- max(tmp.like)
#       log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
#       tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
#     }
#
#     if(out.tau)return(list(logl=log.like,tau=tau))
#     log.like
#   }


# "logLmix_poisson" <- function (pars, first_fit, G){# keep weights and offsets in first fit, weights, offsets, out_tau=FALSE){
#   S <- ncol(first_fit$y)
#   tau <- matrix(0,S,G)
#   if(G>1){
#
#     # remove pi
#     fm <- pars[-1*(1:((G-1)))]  ## remove pi
#
#     # get species intercepts
#     sp_int <- fm[(length(fm)-(S-1)):(length(fm))]
#
#     # get mixture parameters
#     fm <- fm[-1*(length(fm)-(S-1)):(length(fm))]
#
#     # get pi's
#     pi <- pars[(1:(G-1))]
#
#     # re-structure vector to matrix
#     dim(fm) <- c(G,length(fm)/G)
#
#     # re calculate pi's
#     pi <- additive_logistic(pi)
#
#   } else{
#     return(0)
#     fm <- tau[1:(length(pars)-1)]
#     dim(fm) <- c(1,length(fm))
#     pi <- 1
#   }
#
#   log_like <- 0
#   S <- ncol(first_fit$y)
#
#   for(j in 1:S){
#     tmp_like <- rep(0,G)
#     for(i in 1:G){
#       lpre <- first_fit$x%*%c(sp_int[j],fm[i,])+first_fit$offset
#       tmp_like[i] <- sum(dpois(first_fit$y[,j],lambda=exp(lpre),log=TRUE)*first_fit$weights[,j])
#     }
#       eps <- max(tmp_like)
#       log_like <- log_like + (log(sum(pi*exp((tmp_like)-(eps))))+(eps))
#       tau[j,] <- log(pi) + tmp_like - (log(sum(pi*exp((tmp_like)-(eps))))+(eps))
#     }
#   log_like
# }

# "logLmix_tweedie" <-
#   function (pars,first.fit,G,S,sp,spname,out.tau=FALSE)
#   {
#     tau <- matrix(0,S,G)
#     ##tau,out.tau=FALSE
#     if(G>1){
#       pi <- pars[which(names(pars)=="pi")]
#       phi <- pars[which(names(pars)=="phi")]
#       p <- rep(1.6,S)
#       fm <- pars[which(names(pars)=="coef")]
#       sp.int <- pars[which(names(pars)=="int")]
#       ##fm <- pars[-1*(1:((G-1)+G))]  ## remove pi
#       ##sp.int <- fm[(length(fm)-(S*G-1)):length(fm)]
#       ##fm <- fm[-1*((length(fm)-(S*G-1)):length(fm))]
#       ##pi <- pars[(1:(G-1))]
#       ##theta <- pars[G:(G+G-1)]
#
#       ##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
#       #pi <- tau[((length(tau)-(G-2)):length(tau))]
#       ##dim(sp.int) <- c(S,G)
#       dim(fm) <- c(G,length(fm)/G)
#       ##pi[G] <- 1-sum(pi)
#       pi <- additive_logistic(pi)
#
#     } else{
#       return(0)
#       fm <- tau[1:(length(pars)-1)]
#       dim(fm) <- c(1,length(fm))
#       pi <- 1
#     }
#
#
#
#     log.like <- 0
#     for(j in 1:S){
#       sel.sp <- which(sp==spname[j])
#       tmp.like <- rep(0,G)
#       for(i in 1:G){
#         lpre <- cbind(1,first.fit$x[sel.sp,])%*%c(sp.int[j],fm[i,])+first.fit$offset[sel.sp]
#         ##tmp.like[i] <- sum(dpois(first.fit$y[sel.sp],exp(lpre),log=TRUE))
#         tmp.like[i] <- sum(dTweedie(y=first.fit$y[sel.sp],mu=exp(lpre),phi=phi[j],p=p[j],LOG=TRUE))
#       }
#       ##print(tmp.like)
#       eps <- max(tmp.like)
#       log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
#       tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
#
#     }
#     if(out.tau)return(list(logl=log.like,tau=tau))
#     log.like
#   }


# "mix.residuals" <- function (fmM,form,datsp,sp){
#     cat("calculating residuals \n")
#     link.fun <- stats::make.link("logit")
#     x <- stats::model.matrix(form,data=datsp)
#     spname <- unique(sp)
#     S <- length(spname)
#     G <- ncol(fmM$tau)
#
#     PIT <- matrix(NA,S,G)
#     for(g in 1:G){
#       for(s in seq_along(spname)){
#         sel.sp <- which(sp==spname[s])
#         t.obs <- sum(datsp$obs[sel.sp])
#         pre <- link.fun$linkinv(x[sel.sp,]%*%fmM$coef[g,])
#         obs <- datsp$obs[sel.sp]
#         dis <- distr.binom(pre)
#         #      PIT[s,g] <- qnorm(sum(dis[1:(length(dis)-1)],dis[length(dis)]/2))
#         nSucc <- sum( obs)
#         transfo <- sum( dis[1:nSucc],dis[nSucc+1]/2)
#         transfo <- min( transfo, 1)
#         ##transfor <- max( transfo, 0)
#         PIT[s,g] <- qnorm( transfo)
#       }
#     }
#     PIT
#   }

"print.species_mix" <-  function (x,...){
    cat("\nMixing probabilities\n")
    print(x$pi)
    cat("\nCoefficents\n")
    print(x$coef)
    if(!is.na(x$se[1])){
      cat("\nStandard Errors of coefficents\n")
      print(x$se)
    }
    cat("\nPosterior Probabilities\n")
    print(x$tau)
  }


"rPoisGam" <- function ( n, lambda, mu.Z, alpha) {
    mu.N <- lambda
    #simulate n random variables from the same compound poisson distribution
    my.fun <- function (parms)
      return( rgamma( n=1, scale=parms[3], shape=parms[1]*parms[2]))
    #    return( sum( rgamma( n=parms[1], scale=parms[3], shape=parms[2])))

    tmp <- c( n, length( mu.N), length(mu.Z), length( alpha))
    names( tmp) <- c( "n","mu.N","mu.Z","alpha")
    if( !all( is.element( tmp[-1], c( 1, tmp[1])))) {
      print( "rPoisGam: error -- length of arguments are not compatible")
      return( tmp)
    }
    if( tmp["mu.N"]==1)
      mu.N <- rep( mu.N, tmp["n"])
    if( tmp["mu.Z"]==1)
      mu.Z <- rep( mu.Z, tmp["n"])
    if( tmp["alpha"]==1)
      alpha <- rep( alpha, tmp["n"])

    np <- matrix( rpois( n, mu.N), ncol=1)
    beta <- mu.Z / alpha
    y <- apply( cbind( np, alpha, beta), 1, my.fun)

    return( y)
  }

# "species_mix_bernoulli" <- function (sp.form, sp.data,covar.data,G=2, pars=NA, control) {
#     t.covar.data <- covar.data
#     t.sp.data <- sp.data
#     sp.form <- update.formula(sp.form,obs~1+.)
#     if(control$em_prefit | G==1){
#       prefit <- species_mix_em(sp.form,sp.data,covar.data,G,control)
#       pars <- c(additive_logistic(prefit$pi,TRUE)[1:(G-1)],unlist(prefit$coef))
#     }
#     S <- dim(sp.data)[2]
#     if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
#     n <- dim(sp.data)[1]
#     sp <- rep(sp.name,each=n)
#     if(G==1) return(prefit)
#     var.names <- colnames(covar.data)
#     covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
#     names(covar.data) <- var.names
#     sp.data <- as.data.frame(sp.data)
#     data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
#     names(data)[1] <- as.character(sp.form)[2]
#     fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,control)
#     while(fmM.out$logl==0) {
#       prefit <- species_mix_em(sp.form,t.sp.data,t.covar.data,G,control)
#       pars <- c(additive_logistic(prefit$pi,TRUE)[1:(G-1)],unlist(prefit$coef))
#       fmM.out <- fitmix.cpp(sp.form,data,sp,G,pars=pars,control)
#     }
#
#     rownames(fmM.out$tau) <- sp.name ## add names to taus
#     fmM.out$se <- NA
#     if(control$calculate_hessian_cpp){
#       fmM.out$covar <- try(solve(fmM.out$hessian))
#       if(class(fmM.out$covar)!="try-error"){
#         tmp <- sqrt(diag(fmM.out$covar))
#         tmp <- tmp[(G):length(tmp)]
#         fmM.out$se <- matrix(tmp,G,ncol(fmM.out$coef))
#         colnames(fmM.out$se) <- colnames(fmM.out$coef)
#         rownames(fmM.out$se) <- rownames(fmM.out$coef)
#       }
#     }
#     if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
#     fmM.out$formula <- sp.form
#     class(fmM.out) <- c("species_mix","bernoulli")
#     fmM.out
#   }

#
# "species_mix_em" <-  function (sp.form, sp.data, covar.data, G=2, control, residuals=FALSE){
#     S <- dim(sp.data)[2]
#     if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
#     n <- dim(sp.data)[1]
#     sp <- rep(sp.name,each=n)
#
#     var.names <- colnames(covar.data)
#     covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
#     names(covar.data) <- var.names
#     sp.data <- as.data.frame(sp.data)
#     data <- data.frame(obs=unlist(sp.data),covar.data)
#
#     fmM.out <- fitmix(sp.form,data,sp,G,control)
#     if(control$em_refit>1)
#       for(i in 2:control$em_refit){
#         fmM <- fitmix(sp.form,data,sp,G,control)
#         if(fmM$logl>fmM.out$logl) fmM.out <- fmM
#       }
#     if(control$em_calculate_hessian){
#       var <- 0
#       t.pi <- additive_logistic(fmM.out$pi,TRUE)
#       parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef))
#       first.fit <- list(y=data[,1],x=stats::model.matrix(sp.form,data=data))
#       fun_est_var <- function(x){-logLmix(x,first.fit,G,S,sp,sp.name)}
#       var <- solve(numDeriv::hessian(fun_est_var, parms))
#       colnames( var) <- rownames( var) <- names( parms)
#       fmM.out$covar <- var
#     }
#     if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
#     fmM.out$formula <- sp.form
#
#     fmM.out
#   }

##### Bernoulli functions #####

"apply_glmnet_bernoulli_sp" <- function(ss, y, X, weights, offset){
  options(warn=-1)
  lambda.seq <- sort( unique( c(seq( from=1/0.1, to=1, length=10),seq( from=0.1, to=1, length=10))), decreasing=TRUE)
  tmp.fm <- glmnet::glmnet(y=y[,ss], x=X[,-1], weights=weights,
                           offset=offset,
                           family='binomial',
                           alpha=0, #ridge penalty
                           lambda=lambda.seq, #the range of penalties, note that only one will be used
                           standardize=FALSE,  #don't standardize the covariates (they are already standardised)
                           intercept=TRUE)
  my.coefs <- apply(glmnet::coef.glmnet(tmp.fm),1,lambda_penalisation_fun, lambda.seq)
  return(my.coefs)
}

"apply_glmnet_bernoulli_sp_tau" <- function (ss, y, X, offset, tau, G, S, fits){
  Y_sp_tau <- as.matrix(rep(y[,ss],G))
  X_sp_tau <- do.call(rbind, replicate(G, X, simplify=FALSE))
  wts_sp_tau <- rep(tau[ss,],each=length(y[,ss]))
  offy <- rep(offset,G)
  lambda.seq <- sort( unique( c(seq( from=1/0.1, to=1, length=10),seq( from=0.1, to=1, length=10))), decreasing=TRUE)
  f_mix <- glmnet::glmnet(x = X_sp_tau[,-1], y = Y_sp_tau, weights = wts_sp_tau,
                          offset=offy, family='binomial',
                          alpha=0, #ridge penalty
                          lambda=lambda.seq, #the range of penalties, note that only one will be used
                          standardize=FALSE,  #don't standardize the covariates (they are already standardised)
                          intercept=TRUE)
  my.coefs <- apply(glmnet::coef.glmnet(f_mix),1,lambda_penalisation_fun,lambda.seq)
  return(sp_coef=my.coefs)
}

"apply_glmnet_bernoulli_group_tau" <- function (gg, y, X, tau){

  ### setup the data stucture for this model.
  Y_tau <- as.matrix(unlist(as.data.frame(y)))
  X_no_NA <- list()
  for (jj in 1:ncol(y)){
    X_no_NA[[jj]] <- X
  }
  X_tau <- do.call(rbind, X_no_NA)
  n_ys <- sapply(X_no_NA,nrow)
  wts_tau <- rep(tau[,gg],n_ys)

  # ft_mix <- stats::glm.fit(x = X_tau, y = Y_tau, weights = wts_tau, family=stats::binomial())
  lambda.seq <- sort( unique( c(seq( from=1/0.1, to=1, length=10),seq( from=0.1, to=1, length=10))), decreasing=TRUE)
  ft_mix <- glmnet::glmnet(y=Y_tau, x=X_tau[,-1], weights=wts_tau,
                 # offset=offset,
                 family='binomial',
                 alpha=0, #ridge penalty
                 lambda=lambda.seq, #the range of penalties, note that only one will be used
                 standardize=FALSE,  #don't standardize the covariates (they are already standardised)
                 intercept=TRUE)
  mix_coefs <- apply(glmnet::coef.glmnet(ft_mix),1,lambda_penalisation_fun,lambda.seq)[-1]
  # mix_coefs <- ft_mix$coef[-1]
  return(mix_coefs)
}

"get_starting_values_bernoulli_sp" <- function(y, X, weights, offset, G, S, control){
  temp.warn <- getOption( "warn")
  options( warn=-1)

  # S <- ncol(y)
  #if emfit is in the control do an EM fit to get good starting values for c++
  if(isTRUE(control$em_prefit)){
    if(!control$quiet)message('Using EM algorithm to find starting values; using',
                              control$em_refit,'refits\n')
    emfits <- list()
    for(ii in seq_len(control$em_refit)){
      emfits[[ii]] <-  fitmix_EM_bernoulli_sp(y, X, weights, offset, G, control)
    }
    bf <- which.max(vapply(emfits,function(x)c(x$logl),c(logl=0)))
    emfit <- emfits[[bf]]
    start_vals <- list(alphas=emfit$alpha,
                       betas=emfit$beta,
                       pis=emfit$pis)
  } else {
    if(!control$quiet)message('You are not using the EM algorith to find starting values; starting values are
                              generated using',control$init_method,'\n')
    starting_values <- get_initial_values_bernoulli_sp(y = y, X = X, weights = weights,
                                                       offset = offset, G = G, S = S,
                                                       control = control)
    start_vals <- list(alphas=starting_values$fits$alphas,
                       betas=starting_values$fits$betas,
                       pis=starting_values$pis)
  }

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pis,inv = TRUE)[-G]
  start_vals$disp <- rep(NA,S)
  start_vals$spp_wts <- weights
  start_vals$site_spp_wts <- matrix(1,nrow(y),ncol(y))
  start_vals$y_is_na <- matrix(FALSE,nrow(y),ncol(y))
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)

  return(start_vals)
}


"get_logls_bernoulli_sp" <- function(first_fit, fits, G, S){
  logl_sp_bernoulli_sp <- matrix(NA, nrow=S, ncol=G)
  p <- stats::make.link(link = "logit")
  for(ss in 1:S){
    for(gg in 1:G){
      lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
      logl_sp_bernoulli_sp[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, p$linkinv(lp),log = TRUE))
    }
    logl_sp_bernoulli_sp[ss,gg] <- logl_sp_bernoulli_sp[ss,gg]*first_fit$weights[ss]
  }
  return(logl_sp_bernoulli_sp)
}

"get_initial_values_bernoulli_sp" <- function(y, X, weights, offset, G, S, control){#cores, inits='kmeans', init.sd=1){
  starting_values <- initiate_fit_bernoulli_sp(y, X, weights, offset, G, S, control)# cores, inits, init.sd)
  fits <- list(alphas=starting_values$sp_intercepts,betas=starting_values$mix_coefs,disp=rep(NA,S))
  first_fit <- list(x = X, y = y, weights=weights, offset=offset)
  logls <- get_logls_bernoulli_sp(first_fit, fits, G, S)
  pis <- rep(1/G, G)
  taus <- get_taus(pis, logls, G, S)
  taus <- skrink_taus(taus,max_tau=1/G + 0.1, G)

  #use glm for the step.
  fmix_coefs <- surveillance::plapply(1:G, apply_glm_poisson_group_tau,
                                      first_fit$y,
                                      first_fit$x,
                                      taus,
                                      .parallel = control$cores,
                                      .verbose = !control$quiet)

  #update the mix coefs.
  fmix_coefs <- do.call(rbind,fmix_coefs)
  fits$betas <- as.matrix(fmix_coefs)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$pis <- pis
  res$taus <- taus
  return(res)
}


"get_incomplete_logl_bernoulli_sp_function" <-  function(eta, first_fit, fits, G, S){

  p <- stats::make.link('logit')
  pis <- additive_logistic(eta)
  logl_sp_ippm <- matrix(NA,nrow=S, ncol=G)
  for(ss in 1:S){
    for(gg in 1:G){
      #lp is the same as log_mu (linear predictor)
      lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
      logl_sp_ippm[ss,gg] <- sum(dbinom(first_fit$y[,ss], 1, p$linkinv(lp),log = TRUE)*first_fit$weights)
    }
  }
  # print(log(pis))
  ak <- logl_sp_ippm + matrix(rep(log(pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  logl <- sum( sppLogls)
  return( logl)
}


## gradient stuff
# "incom_logl_bernoulli_sp_alpha" <- function(x, first_fit, eta, fits, G, S){
#   fits$alphas <- x
#   # pis <- additive_logistic(eta)
#   tmp <- get_incomplete_logl_bernoulli_sp_function(eta, first_fit, fits, G, S)
#   return(-tmp)
# }

"incom_logl_bernoulli_sp_beta" <- function(x, first_fit, eta, fits, G, S){

  fits$betas <- matrix(x,nrow=nrow(fits$betas),ncol=ncol(fits$betas))
  # pis <- additive_logistic(eta)
  tmp <- get_incomplete_logl_bernoulli_sp_function(eta, first_fit, fits, G, S)
  return(-tmp)
}

# "incom_logl_bernoulli_sp_pi" <- function(x, first_fit, fits, G, S){
#   eta <- x
#   # pis <- additive_logistic(eta)
#   tmp <- get_incomplete_logl_bernoulli_sp_function(eta, first_fit, fits, G, S)
#   return(-tmp)
# }

"initiate_fit_bernoulli_sp" <- function(y, X, weights, offset, G, S, control){#cores, inits='kmeans', init.sd=1){
  fm_bern_sp_int <- surveillance::plapply(1:S, apply_glmnet_bernoulli_sp, y, X, weights, offset, .parallel = control$cores, .verbose = !control$quiet)
  all_coefs <- do.call(rbind, fm_bern_sp_int)
  mix_coefs <- all_coefs[,-1] # drop intercepts
  # print(all_coefs)
  if(control$init_method=='kmeans'){
    if(!control$quiet)cat("Initial groups by K-means clustering\n")
    tmp1 <- stats::kmeans( mix_coefs, centers=G, nstart=100)
    tmp_grp <- tmp1$cluster
    grp_coefs <- apply( mix_coefs, 2, function(x) tapply(x, tmp_grp, mean))
  }
  if(control$init_method=='hclust'){
    if(!control$quiet)cat("Initial groups by hierarchical clustering using Ward's method on Euclidean distances\n")
    disto <- dist( mix_coefs, method="euclidean")
    tmp1 <- hclust(disto, method = "ward.D")
    tmp_grp <- cutree(tmp1, G)
    grp_coefs <- apply( mix_coefs, 2, function(x) tapply(x, tmp_grp, mean))
  }
  if(control$init_method=='random'){
    if(!control$quiet)cat("Initial groups by random allocation and means from random numbers\n")
    grp_coefs <- matrix( stats::rnorm(G*ncol(mix_coefs), sd=control$init_sd, mean=0), nrow=G, ncol=ncol(mix_coefs))
    tmp_grp <- sample(1:G, S, replace=TRUE)
  }
  colnames(grp_coefs)<-colnames(X[,-1])
  results <- list()
  results$grps <- tmp_grp
  results$mix_coefs <- grp_coefs
  results$sp_intercepts <- all_coefs[,1]
  results$all_coefs <- all_coefs
  return(results)
}

"fitmix_EM_bernoulli_sp" <- function(y, X, weights, offset, G, control){

  S <- ncol(y)
  pis <- rep(0, G)
  ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888

  # get starting values
  starting_values <- get_initial_values_bernoulli_sp(y = y, X = X,
                                                     weights = weights,
                                                     offset = offset,
                                                     G = G, S = S,
                                                     control = control)

  # first e-step
  pis <- starting_values$pis
  fits <- starting_values$fits
  taus <- starting_values$taus
  first_fit <- starting_values$first_fit

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){

    # Estimate pis from tau.
    pis <- colSums(taus)/S

    if (any(pis == 0)) {
      starting_values <- get_initial_values_bernoulli_sp(y = y, X = X,
                                                         weights = weights, offset = offset,
                                                         G = G, S = S,
                                                         control = control)

      pis <- starting_values$pis
      fits <- starting_values$fits
      taus <- starting_values$taus
      first_fit <- starting_values$first_fit
      ite <- 1
    }
    # m-step
    tmp <- stats::nlminb(start=fits$betas, objective=incom_logl_bernoulli_sp_beta, gradient=NULL, hessian=NULL,
                         first_fit=first_fit, eta=additive_logistic(pis,inv = TRUE)[-G], fits=fits, G=G, S=S)
    fits$betas <- update_mix_coefs(fits$betas, tmp$par)
    fm_bernoulli_sp_int <- surveillance::plapply(1:S, apply_glmnet_bernoulli_sp, y, X, weights, offset, .parallel = control$cores, .verbose = !control$quiet) #check weights in this.
    sp_int <- do.call(cbind,fm_bernoulli_sp_int)[1,]
    fits$alphas <- update_sp_coefs(fits$alphas,sp_int)

    # e-step
    # get the log-likes and taus
    logls <- get_logls_bernoulli_sp(first_fit, fits, G, S)
    taus <- get_taus(pis, logls, G, S)

    #update the likelihood
    logl_old <- logl_new
    logl_new <- get_incomplete_logl_bernoulli_sp_function(eta = additive_logistic(pis,inv = TRUE)[-G], first_fit, fits, G, S)
    ite <- ite + 1
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", 1:G, sep = "")
  int_out <- fits$alphas
  fm_out <- fits$betas
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- additive_logistic(pis, TRUE)[-1]

  # estimate log-likelihood
  logl_new <- get_incomplete_logl_bernoulli_sp_function(eta, first_fit, fits, G, S)

  return(list(logl = logl_new, alpha = int_out, beta = fm_out,
              eta = eta, pis = pis, taus = round(taus,4)))

}

###### Poisson functions ####

# "apply_glm_poisson_sp_tau" <- function (ss, y, X, offset, tau, G, S, fits){
#   Y_sp_tau <- as.matrix(rep(y[,ss],G))
#   X_sp_tau <- do.call(rbind, replicate(G, X, simplify=FALSE))
#   wts_sp_tau <- rep(tau[ss,],each=length(y[,ss]))
#   offy <- rep(offset,G)
#   lambda.seq <- sort( unique( c(seq( from=1/0.1, to=1, length=10),seq( from=0.1, to=1, length=10))), decreasing=TRUE)
#   f_mix <- glmnet::glmnet(x = X_sp_tau[,-1], y = Y_sp_tau, weights = wts_sp_tau,
#                           offset=offy, family='poisson',
#                           alpha=0, #ridge penalty
#                           lambda=lambda.seq, #the range of penalties, note that only one will be used
#                           standardize=FALSE,  #don't standardize the covariates (they are already standardised)
#                           intercept=TRUE)
#   my.coefs <- apply(glmnet::coef.glmnet(f_mix),1,lambda_penalisation_fun)
#   return(sp_coef=my.coefs)
# }
#

"apply_glmnet_poisson" <- function(ss, y, X, weights, offset){

  lambda.seq <- sort( unique( c( seq( from=1/0.01, to=1, length=10), seq( from=1/0.1, to=1, length=10),seq(from=0.9, to=10^-2, length=10))), decreasing=TRUE)
  # ids_i <- !y_is_na[,i]
  ft_pos <- glmnet::glmnet(x=as.matrix(X[,-1]),y=as.matrix(y[,ss]),weights=weights,
                           offset=offset,family='poisson',
                           alpha=0,
                           lambda = lambda.seq,
                           standardize = FALSE,
                           intercept = TRUE)
  locat_s <- min(lambda.seq) # relatively small penalty 1/10 - the other ones were to restrictive.
  my_coefs <- glmnet::coef.glmnet(ft_pos, s=locat_s)
  if( any( is.na( my_coefs))){  #just in case the model is so badly posed that mild penalisation doesn't work...
    my_coefs <- glmnet::coef.glmnet(ft_poisson, s=lambda.seq)
    lastID <- apply( my_coefs, 2, function(x) !any( is.na(x)))
    lastID <- tail( (seq_along( lastID))[lastID], 1)
    my_coefs <- my_coefs[,lastID]
  }
  return(as.matrix(my_coefs))
}

"apply_glm_poisson_group_tau" <- function (gg, y, X, tau){

  ### setup the data stucture for this model.
  Y_tau <- as.matrix(unlist(as.data.frame(y)))
  X_no_NA <- list()
  for (jj in 1:ncol(y)){
    X_no_NA[[jj]] <- X
  }
  X_tau <- do.call(rbind, X_no_NA)
  n_ys <- sapply(X_no_NA,nrow)
  wts_tau <- rep(tau[,gg],n_ys)

  ft_mix <- stats::glm.fit(x = X_tau, y = Y_tau, weights = wts_tau, family=stats::poisson())
  mix_coefs <- ft_mix$coef[-1]
  return(mix_coefs)
}

"get_starting_values_poisson" <- function(y, X, weights, offset, G, S, control){ #try and keep control at the end

  if(isTRUE(control$em_prefit)){
    if(!control$quiet)message('Using EM algorithm to find starting values; using ',
                              control$em_refit,' refits\n')
    emfits <- list()
    for(ii in seq_len(control$em_refit)){
      emfits[[ii]] <-  fitmix_EM_poisson(y, X, weights, offset, G, S, control)
    }
    bf <- which.max(vapply(emfits,function(x)c(x$logl),c(logl=0)))
    emfit <- emfits[[bf]]
    start_vals <- list(alphas=emfit$alpha,
                       betas=emfit$beta,
                       pis=emfit$pis)
  } else {
    if(!control$quiet)message('You are not using the EM algorith to find starting values; starting values are
                              generated using',control$init_method,'\n')

    starting_values <- get_initial_values_poisson(y = y, X = X, weights = weights,
                                               offset = offset,
                                               G = G, S = S,
                                               control = control)
    start_vals <- list(alphas=starting_values$fits$alphas,
                       betas=starting_values$fits$betas,
                       pis=starting_values$pis)
  }

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pis,inv = TRUE)[-G]
  start_vals$disp <- rep(NA,S)
  start_vals$spp_wts <- weights # bootstrap weights - not doing anything here.
  start_vals$site_spp_wts <- matrix(1,nrow(y),ncol(y)) #ippm weights
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)
  return(start_vals)
}

## new version with glmnet for regularisation.
"get_initial_values_poisson" <- function(y, X, weights, offset, G, S, control){

  # fit glmnet to get the intial values
  starting_values <- initiate_fit_poisson(y, X, weights, offset, G, S, control)
  fits <- list(betas=starting_values$mix_coefs, alphas=starting_values$sp_intercepts)
  first_fit <- list(x = X, y = y, weights=weights, offset=offset)

  # get the loglikelihood based on these values
  logls <- get_logls_poisson(first_fit, fits, G, S)

  # estimate the posteriors for taus (skrink a little)
  pis <- rep(1/G, G)
  taus <- get_taus(pis, logls, G, S)
  taus <- skrink_taus(taus, max_tau=1/G + 0.1, G)

  #use glm for the step.
  fmix_coefs <- surveillance::plapply(1:G, apply_glm_poisson_group_tau,
                                      first_fit$y,
                                      first_fit$x,
                                      taus,
                                      .parallel = control$cores,
                                      .verbose = !control$quiet)

  #update the mix coefs.
  fmix_coefs <- do.call(rbind,fmix_coefs)
  fits$betas <- as.matrix(fmix_coefs)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$pis <- pis
  res$taus <- taus
  return(res)
}

"get_incomplete_logl_poisson" <-  function(eta, first_fit, fits, G, S){
  pis <- additive_logistic(eta)
  logl_sp_poisson <- matrix(NA, nrow=S, ncol=G)
  for(ss in 1:S){
    for(gg in 1:G){
      #eta is the same as log_lambda (linear predictor)
      eta <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
      logl_sp_poisson[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(eta),log=TRUE))
    }
  }
  ak <- logl_sp_poisson + matrix( rep( log( pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  logl <- sum( sppLogls)
  return( logl)
}

"get_logls_poisson" <- function(first_fit, fits, G, S){
  logl_sp_poisson <- matrix(NA, nrow=S, ncol=G)
  for(ss in 1:S){
    # sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      #lp is the same as log_lambda (linear predictor)
      lp <- first_fit$x[,1] * fits$alphas[ss] + as.matrix(first_fit$x[,-1]) %*% fits$betas[gg,] + first_fit$offset
      logl_sp_poisson[ss,gg] <- sum(dpois(first_fit$y[,ss],exp(lp),log=TRUE))
    }
  }
  return(logl_sp_poisson)
}


"initiate_fit_poisson" <- function(y, X, weights, offset, G, S, control){
  fm_pos <- surveillance::plapply(1:S,apply_glmnet_poisson, y, X, weights, offset,
                                   .parallel = control$cores, .verbose = !control$quiet)
  # fm_pos <- surveillance::plapply(1:S,apply_glm_poisson, y, X, weights, offset, y_is_na,
  #                                  .parallel = control$cores, .verbose = !control$quiet)
  all_coefs <- t(do.call(cbind,fm_pos)) #be careful with this transformation from glmnet.
  mix_coefs <- all_coefs[,-1] # drop intercepts
  if(control$init_method=='kmeans'){
    if(!control$quiet)message( "Initial groups by K-means clustering\n")
    tmp1 <- stats::kmeans(mix_coefs, centers=G, nstart=100)
    tmp_grp <- tmp1$cluster
    grp_coefs <- apply(mix_coefs, 2, function(x) tapply(x, tmp_grp, mean))
  }
  if(control$init_method=='hclust'){
    stop( "EHNT!... Poisson cannot use a hierachial cluster - stick to K-means\n")
  }
  if(control$init_method=='random' | is.null(tmp_grp)){
    if(!control$quiet)message( "Initial groups by random allocation and means from random numbers\n")
    grp_coefs <- matrix( stats::rnorm(G*ncol(mix_coefs), sd=control$init_sd, mean=0), nrow=G, ncol=ncol(mix_coefs))
    tmp_grp <- sample(1:G, S, replace=TRUE)
  }

  colnames(grp_coefs) <- colnames(X[,-1])
  results <- list()
  results$grps <- tmp_grp
  results$mix_coefs <- grp_coefs
  results$sp_intercepts <- all_coefs[,1]
  results$all_coefs <- all_coefs
  return(results)
}

"incom_logl_poisson" <- function(x, first_fit, pis, fits, G, S){
  fits$betas <- matrix(x,nrow=nrow(fits$betas),ncol=ncol(fits$betas))
  tmp <- get_incomplete_logl_poisson(pis, first_fit, fits, G, S)
  return(-tmp)
}

"incom_logl_poisson_beta" <- function(x, first_fit, eta, fits, G, S){

  fits$betas <- matrix(x,nrow=nrow(fits$betas),ncol=ncol(fits$betas))
  # pis <- additive_logistic(eta)
  tmp <- get_incomplete_logl_poisson(eta, first_fit, fits, G, S)
  return(-tmp)
}

"fitmix_EM_poisson" <- function(y, X, weights, offset, G, S, control){

  # S <- ncol(y)
  pis <- rep(0, G)
  ite <- 1
  logl_old <- -99999999
  logl_new <- -88888888
  control$quiet <- TRUE
  # get starting values
  starting_values <- get_initial_values_poisson(y = y, X = X, offset = offset,
                                                     weights = weights, G = G, S = S,
                                                     control = control)

  # first e-step
  pis <- starting_values$pis
  fits <- starting_values$fits
  taus <- starting_values$taus
  first_fit <- starting_values$first_fit

  while(control$em_reltol(logl_new,logl_old) & ite <= control$em_steps){

    # Estimate pis from tau.
    pis <- colSums(taus)/S

    if (any(pis == 0)) {
      starting_values <- get_initial_values_poisson(y = y, X = X, offset = offset,
                                                         weights = weights, G = G, S = S,
                                                         control = control)

      pis <- starting_values$pis
      fits <- starting_values$fits
      taus <- starting_values$taus
      first_fit <- starting_values$first_fit
      ite <- 1
    }
    # m-step
    tmp <- stats::nlminb(start=fits$betas, objective=incom_logl_poisson_beta, gradient=NULL, hessian=NULL,
                         first_fit=first_fit, eta=additive_logistic(pis,inv = TRUE)[-G], fits=fits, G=G, S=S)
    fits$betas <- update_mix_coefs(fits$betas, tmp$par)
    fm_poissonint <- surveillance::plapply(1:S, apply_glmnet_poisson, y, X, weights, offset, .parallel = control$cores, .verbose = !control$quiet) #check weights in this.
    sp_int <- do.call(cbind,fm_poissonint)[1,]
    fits$alphas <- update_sp_coefs(fits$alphas,sp_int)

    # e-step
    # get the log-likes and taus
    logls <- get_logls_poisson(first_fit, fits, G, S)
    taus <- get_taus(pis, logls, G, S)

    #update the likelihood
    logl_old <- logl_new
    logl_new <- get_incomplete_logl_poisson(eta = additive_logistic(pis,inv = TRUE)[-G], first_fit, fits, G, S)
    ite <- ite + 1
  }

  taus <- data.frame(taus)
  names(taus) <- paste("grp.", 1:G, sep = "")
  int_out <- fits$alphas
  fm_out <- fits$betas
  names(pis) <- paste("G", 1:G, sep = ".")
  eta <- additive_logistic(pis, TRUE)[-G]

  # estimate log-likelihood
  logl_new <- get_incomplete_logl_poisson(eta, first_fit, fits, G, S)

  return(list(logl = logl_new, alpha = int_out, beta = fm_out,
              eta = eta, pis = pis, taus = round(taus,4)))

}

###### IPPM functions #######


# "apply_glm_ippm" <- function(i, y, X, weights, offset, y_is_na){
#   ids_i <- !y_is_na[,i]
#   f_ippm <- stats::glm.fit(x=X[ids_i,],y=y[ids_i,i]/weights[ids_i,i],weights=weights[ids_i,i],offset=offset[ids_i],family=stats::poisson())
#   f_ippm$coef
# }

"apply_glmnet_ippm" <- function(i, y, X, weights, offset, y_is_na,return_all_coefs=FALSE){

  lambda.seq <- sort( unique( c( seq( from=1/0.01, to=1, length=10), seq( from=1/0.1, to=1, length=10),seq(from=0.9, to=10^-2, length=10))), decreasing=TRUE)
  ids_i <- !y_is_na[,i]
  ft_ippm <- glmnet::glmnet(x=as.matrix(X[ids_i,-1]),y=as.matrix(y[ids_i,i]/weights[ids_i,i]),weights=as.matrix(weights[ids_i,i]),
                            offset=offset[ids_i],family='poisson',
                            alpha=0,
                            lambda = lambda.seq,
                            standardize = FALSE,
                            intercept = TRUE)
  locat_s <- min(lambda.seq) # relatively small penalty 1/10 - the other ones were to restrictive.
  my_coefs <- glmnet::coef.glmnet(ft_ippm, s=locat_s)
  if(return_all_coefs)  my_coefs <- glmnet::coef.glmnet(ft_ippm)
  if( any( is.na( my_coefs))){  #just in case the model is so badly posed that mild penalisation doesn't work...
    my_coefs <- glmnet::coef.glmnet(ft_ippm, s=lambda.seq)
    lastID <- apply( my_coefs, 2, function(x) !any( is.na(x)))
    lastID <- tail( (seq_along( lastID))[lastID], 1)
    my_coefs <- my_coefs[,lastID]
  }
  return(as.matrix(my_coefs))
}

"apply_glm_ippm_group_tau" <- function (gg, y, X, y_is_na, tau){

  ### setup the data stucture for this model.
  Y_tau <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
  X_no_NA <- list()
  for (jj in 1:ncol(y)){
    X_no_NA[[jj]] <- X[!y_is_na[,jj],]
  }
  X_tau <- do.call(rbind, X_no_NA)
  n_ys <- sapply(X_no_NA,nrow)
  wts_tau <- rep(tau[,gg],n_ys)

  ft_mix <- stats::glm.fit(x = X_tau, y = Y_tau, weights = wts_tau, family=stats::poisson())
  mix_coefs <- ft_mix$coef[-1]
  return(mix_coefs)
}

# "apply_glmnet_ippm_group_tau_v2" <- function (gg, y, X, weights, offset, y_is_na, tau, lambda_pen = NULL, return_all_coefs=FALSE){#currently no weights yet.
#
#   ### setup the data stucture for this model.
#   Y_tau <- as.matrix(unlist(as.data.frame(y[!y_is_na])))
#   X_no_NA <- list()
#   for (jj in 1:ncol(y)){
#     X_no_NA[[jj]] <- X[!y_is_na[,jj],]
#   }
#   X_tau <- do.call(rbind, X_no_NA)
#   n_ys <- sapply(X_no_NA,nrow)
#   wts_tau <- rep(tau[,gg],c(n_ys))
#
#
#   ippm_weights <- as.matrix(as.matrix(unlist(as.data.frame(weights[!y_is_na]))))
#   Z_tau <- as.matrix(Y_tau/ippm_weights)
#   wts_tauXippm_weights <- wts_tau*ippm_weights
#   offy_mat <- replicate(ncol(y),offset)
#   offy <- unlist(as.data.frame(offy_mat[!y_is_na]))
#
#   ## setup the data penality parameters for glmnet
#   # lambda.seq <- 10^seq(1,-2,length=50)
#   lambda.seq <- sort( unique( c( seq( from=1/0.01, to=1, length=10), seq( from=1/0.1, to=1, length=10),seq(from=0.9, to=10^-2, length=10))), decreasing=TRUE)
#   ft_mix <- glmnet::glmnet(x=as.matrix(X_tau[,-1]),y=as.matrix(Z_tau),
#                            weights=as.matrix(wts_tauXippm_weights),offset = offy,
#                            family='poisson',
#                            alpha=0,
#                            lambda = lambda.seq,
#                            standardize = FALSE,
#                            intercept = FALSE)
#   if(is.null(lambda_pen))lambda_pen <- min(lambda.seq) # relatively small penalty 1/10 - the other ones were to restrictive
#   my_coefs <- glmnet::coef.glmnet(ft_mix, s=lambda_pen)
#   if(return_all_coefs)  my_coefs <- glmnet::coef.glmnet(ft_mix)
#   if(any(is.na( my_coefs))){  #just in case the model is so badly posed that mild penalisation doesn't work...
#     my_coefs <- glmnet::coef.glmnet(ft_mix, s=lambda.seq)
#     lastID <- apply( my_coefs, 2, function(x) !any( is.na(x)))
#     lastID <- tail( (seq_along( lastID))[lastID], 1)
#     my_coefs <- my_coefs[,lastID]
#   }
#   return(as.matrix(my_coefs))
# }

"get_starting_values_ippm" <- function(y, X, weights, offset, S, G, y_is_na, control){ #try and keep control at the end

  starting_values <- get_initial_values_ippm(y = y, X = X, weights = weights,
                                             offset = offset, y_is_na = y_is_na,
                                             G = G, S = S,
                                             control = control)
  start_vals <- list(alphas=starting_values$fits$alphas,
                     betas=starting_values$fits$betas,
                     pis=starting_values$pis)

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pis,inv = TRUE)[-G]
  start_vals$disp <- rep(NA,S)
  start_vals$spp_wts <- rep(1,nrow(y)) # bootstrap weights - not doing anything here.
  start_vals$site_spp_wts <- weights #ippm weights
  start_vals$y_is_na <- y_is_na #which data are missing.
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)
  return(start_vals)
}


## new version with glmnet for regularisation.
"get_initial_values_ippm" <- function(y, X, weights, offset, y_is_na, G, S, control){

  # fit glmnet to get the intial values
  starting_values <- initiate_fit_ippm(y, X, weights, offset, y_is_na, G, S, control)
  fits <- list(betas=starting_values$mix_coefs, alphas=starting_values$sp_intercepts)
  first_fit <- list(x = X, y = y, weights=weights, offset=offset, y_is_na=y_is_na)

  # get the loglikelihood based on these values
  logls <- get_logls_ippm(first_fit, fits, G, S)

  # estimate the posteriors for taus (skrink a little)
  pis <- rep(1/G, G)
  taus <- get_taus(pis, logls, G, S)
  # taus <- skrink_taus(taus, max_tau=1/G + 0.1, G)

  # use these posteriors to estimate mix-coefs again with weights
  # could replace with glmnet if I can get it to work.
  # apply_glm_ippm_group_tau()
  # fmix_coefs <- surveillance::plapply(1:G, apply_glmnet_ippm_group_tau_v2,
  #                                     first_fit$y,
  #                                     first_fit$x,
  #                                     first_fit$weights,
  #                                     first_fit$offset,
  #                                     first_fit$y_is_na,
  #                                     taus,
  #                                     .parallel = control$cores,
  #                                     .verbose = !control$quiet)
  #
  # #update the mix coefs.
  # fmix_coefs <- t(do.call(cbind,fmix_coefs))
  # fits$betas <- as.matrix(fmix_coefs[,-1])

  #use glm for the step.
  fmix_coefs <- surveillance::plapply(1:G, apply_glm_ippm_group_tau,
                                      first_fit$y,
                                      first_fit$x,
                                      # first_fit$weights,
                                      # first_fit$offset,
                                      first_fit$y_is_na,
                                      taus,
                                      .parallel = control$cores,
                                      .verbose = !control$quiet)

  #update the mix coefs.
  fmix_coefs <- do.call(rbind,fmix_coefs)
  fits$betas <- as.matrix(fmix_coefs)

  res <- list()
  res$fits <- fits
  res$first_fit <- first_fit
  res$pis <- pis
  res$taus <- taus
  return(res)
}

"get_incomplete_logl_ippm" <-  function(pis, first_fit, fits, G, S){
  incomplete.logl <- 0
  logl_sp_ippm <- matrix(NA, nrow=S, ncol=G)
  for(ss in 1:S){
    sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      #eta is the same as log_lambda (linear predictor)
      lp <- first_fit$x[sp_idx,1] * fits$alphas[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$betas[gg,] + first_fit$offset[sp_idx]
      logl_sp_ippm[ss,gg] <- first_fit$y[sp_idx,ss] %*% lp - first_fit$weights[sp_idx,ss] %*% exp(lp)
    }
  }
  ak <- logl_sp_ippm + matrix( rep( log( pis), each=S), nrow=S, ncol=G)
  am <- apply( ak, 1, max)
  ak <- exp( ak-am)
  sppLogls <- am + log( rowSums( ak))
  logl <- sum( sppLogls)
  return( logl)
}

"get_logls_ippm" <- function(first_fit, fits, G, S){
  logl_sp_ippm <- matrix(NA, nrow=S, ncol=G)
  for(ss in 1:S){
    sp_idx<-!first_fit$y_is_na[,ss]
    for(gg in 1:G){
      #lp is the same as log_lambda (linear predictor)
      lp <- first_fit$x[sp_idx,1] * fits$alphas[ss] + as.matrix(first_fit$x[sp_idx,-1]) %*% fits$betas[gg,] + first_fit$offset[sp_idx]
      logl_sp_ippm[ss,gg] <- first_fit$y[sp_idx,ss] %*% lp - first_fit$weights[sp_idx,ss] %*% exp(lp)
    }
  }
  return(logl_sp_ippm)
}

"incom_logl_ippm" <- function(x, first_fit, pis, fits, G, S){
  fits$betas <- matrix(x,nrow=nrow(fits$betas),ncol=ncol(fits$betas))
  tmp <- get_incomplete_logl_ippm(pis, first_fit, fits, G, S)
  return(-tmp)
}

"initiate_fit_ippm" <- function(y, X, weights, offset, y_is_na, G, S, control){
  fm_ippm <- surveillance::plapply(1:S,apply_glmnet_ippm, y, X, weights, offset, y_is_na,
                                   .parallel = control$cores, .verbose = !control$quiet)
  # fm_ippm <- surveillance::plapply(1:S,apply_glm_ippm, y, X, weights, offset, y_is_na,
  #                                  .parallel = control$cores, .verbose = !control$quiet)
  all_coefs <- t(do.call(cbind,fm_ippm)) #be careful with this transformation from glmnet.
  mix_coefs <- all_coefs[,-1] # drop intercepts
  if(control$init_method=='kmeans'){
    if(!control$quiet)message( "Initial groups by K-means clustering\n")
    tmp1 <- stats::kmeans(mix_coefs, centers=G, nstart=100)
    tmp_grp <- tmp1$cluster
    grp_coefs <- apply(mix_coefs, 2, function(x) tapply(x, tmp_grp, mean))
  }
  if(control$init_method=='hclust'){
    stop( "EHNT!... IPPM cannot use a hierachial cluster - stick to K-means\n")
  }
  if(control$init_method=='random' | is.null(tmp_grp)){
    if(!control$quiet)message( "Initial groups by random allocation and means from random numbers\n")
    grp_coefs <- matrix( stats::rnorm(G*ncol(mix_coefs), sd=control$init_sd, mean=0), nrow=G, ncol=ncol(mix_coefs))
    tmp_grp <- sample(1:G, S, replace=TRUE)
  }

  colnames(grp_coefs) <- colnames(X[,-1])
  results <- list()
  results$grps <- tmp_grp
  results$mix_coefs <- grp_coefs
  results$sp_intercepts <- all_coefs[,1]
  results$all_coefs <- all_coefs
  return(results)
}



##### Negative Binomial Functions #####
"apply_glm_nbinom" <- function (i,form,datsp,tau,n){
  dat.tau <- rep(tau[,i],each=n)
  x <- stats::model.matrix(stats::as.formula(form),data=datsp)
  y <- datsp$obs
  environment(form) <- environment()
  f.mix <- nbglm(stats::as.formula(form),data=datsp,weights=dat.tau)
  sp.int <- rep(f.mix$coef[1],dim(tau)[1])
  return(list(coef=f.mix$coef[-1],theta=f.mix$theta,sp.intercept=sp.int))
}

"create_starting_values_nbinom" <- function (S,G,n,form,datsp,control) {
  environment(form) <- environment()
  tau <- matrix(stats::runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
    pi[i] <- sum(tau[,i])/S
  }
  fmM <- surveillance::plapply(1:G,apply_glm_nbinom,form,datsp,tau,n,.parallel=control$cores, .verbose = !control$quiet)
  offset <- stats::model.frame(stats::as.formula(form),data=datsp)
  offset <- stats::model.offset(offset)
  if(is.null(offset)) offset <- rep(0,length(datsp$obs))
  first.fit <- list(x=stats::model.matrix(stats::as.formula(form),data=datsp)[,-1],y=datsp$obs,offset=offset,formula=form)
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}


"create_starting_values_nbinom_kmeans" <- function (S, G, n, form, datsp, tol = 0.1, control){
  MM <- stats::model.matrix(form, datsp)
  offset <- stats::model.frame(form,datsp)
  offset <- stats::model.offset(offset)
  if(is.null(offset)) offset <- rep(0,nrow(MM))

  first.fit <- list(x = stats::model.matrix(stats::as.formula(form), data = datsp)[,
                                                                                   -1], y = datsp$obs, offset=offset, formula = form)
  if(class(first.fit$x)=="numeric") {first.fit$x <- matrix(first.fit$x,length(first.fit$x),1) }#deal with model matrix returning a vector for covar == 1

  if (tol < 0 || tol >= 1)
    stop("Minimum Prevalence % must be between 0 and 1")
  sp.name <- 1:S
  sp <- rep(sp.name, each = n)
  starting.fitem <- list(intercepts = rep(0, S), alpha = rep(0,
                                                             S))
  all.betas <- matrix(0, nrow = S, ncol = ncol(MM))
  colnames(all.betas) <- colnames(MM)
  for (j in 1:S) {
    fit <- nbglm(form, datsp[((j - 1) * n):(j * n - 1), ], est_var = FALSE)
    starting.fitem$sp.intercepts[j] <- fit$coef[1]
    all.betas[j, ] <- fit$coef
    starting.fitem$theta[j] = fit$theta
  }
  all.betas[, 1] <- 0
  cat("Clustering...\n")
  fmmvnorm <- stats::kmeans(x = all.betas, centers = G, iter.max = 100,
                            nstart = 50)
  starting.fitem$coef <- fmmvnorm$centers
  fmM <- list()
  for (i in 1:G) {
    B <- matrix(rep(fmmvnorm$centers[i, ], nrow(datsp)),
                nrow(datsp), ncol(fmmvnorm$centers), byrow = TRUE)
    B[, 1] <- rep(starting.fitem$sp.intercepts, each = n)
    fitted <- exp(rowSums(MM * B)+offset)
    fmM[[i]] <- list(coef = fmmvnorm$centers[i, 2:ncol(fmmvnorm$centers)],
                     theta = mean(starting.fitem$theta[fmmvnorm$cluster ==
                                                         i]), sp.intercept = starting.fitem$sp.intercepts,
                     fitted = fitted)
  }
  tau <- matrix(0, S, G)
  pi <- rep(1/G, G)
  pi <- stats::runif(G, 0.2, 0.8)
  pi <- pi/sum(pi)
  est.tau <- surveillance::plapply(1:S, estimate_pi_nbinom, sp, sp.name, datsp,
                                   fmM, pi, G, first.fit,.parallel=control$cores, .verbose = !control$quiet)
  max.newTau <- 0.8
  alpha <- (1 - max.newTau * G)/(max.newTau * (2 - G) - 1)
  for (j in 1:S) {
    newTau <- (2 * alpha * est.tau[[j]]$tau - alpha + 1)/(2 *
                                                            alpha - alpha * G + G)
    tau[j, ] <- newTau
  }
  for (i in 1:G) {
    pi[i] <- sum(tau[, i])/S
  }
  return(list(pi = pi, fmM = fmM, tau = tau, first.fit = first.fit))
}

"estimate_pi_nbinom" <- function (j,sp,spname,datsp,fmM,pi,G,first.fit){

  tmp.like <- rep(0,G)
  tau <- rep(0,G)
  sel.sp <- which(sp==spname[j])

  for(i in 1:G) {
    lpre <- stats::dnbinom(first.fit$y[sel.sp],mu=fmM[[i]]$fitted[sel.sp],size=fmM[[i]]$theta,log=TRUE)
    tmp.like[i] <- sum(lpre)
  }

  eps <- max(tmp.like)
  sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  for(i in 1:G) {
    tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
  }
  return(list(tau=tau,sum.like=sum.like))
}

"fitmix_nbinom" <- function (form, datsp, sp, G = 2, control){
  temp.warn <- getOption("warn")
  options(warn = -1)
  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp == sp.name[1]))
  cat("Fitting Group", G, "\n")
  if (trace==1)
    cat("Iteration | LogL \n")
  dat.tau <- 0
  pi <- rep(0, G)
  ite <- 1
  logL <- -99999999
  old.logL <- -88888888
  if (ite != 1)
    t1 <- create_starting_values_nbinom(S, G, n, form, datsp)
  t1 <- create_starting_values_nbinom_kmeans(S, G, n, form, datsp)
  pi <- t1$pi
  fmM <- t1$fmM
  tau <- t1$tau
  dat.tau <- t1$dat.tau
  first.fit <- t1$first.fit
  while (abs(logL - old.logL) > 1e-04 & ite <= control$maxit) {
    old.logL <- logL
    for (i in 1:G) {
      pi[i] <- sum(tau[, i])/S
    }
    if (any(pi == 0)) {
      cat("pi has gone to zero - restarting fitting \n")
      t1 <- create_starting_values_nbinom(S, G, n, form, datsp)
      pi <- t1$pi
      fmM <- t1$fmM
      tau <- t1$tau
      dat.tau <- t1$dat.tau
      first.fit <- t1$first.fit
      ite <- 1
    }
    fmM <- surveillance::plapply(1:G, weighted_glm_nbinom, first.fit, tau, n, fmM, sp,.parallel=control$cores, .verbose = !control$quiet)
    for (j in 1:S) {
      tmp <- rep(0, G)
      for (g in 1:G) tmp[g] <- fmM[[g]]$sp.intercept[j]
      tmp <- sum(tmp * tau[j, ])
      for (g in 1:G) fmM[[g]]$sp.intercept[j] <- tmp
    }
    logL <- 0
    tmp.like <- matrix(0, S, G)
    est.tau <- surveillance::plapply(1:S, estimate_pi_nbinom, sp, sp.name, datsp, fmM, pi, G, first.fit, .parallel=control$cores, .verbose = !control$quiet)
    for (j in 1:S) {
      if (is.atomic(est.tau[[j]])) {
        print(est.tau[[j]])
      }
      else {
        tau[j, ] <- est.tau[[j]]$tau
        logL <- logL + est.tau[[j]]$sum.like
      }
    }
    if (trace)
      cat(ite, "     | ", logL, "\n")
    ite <- ite + 1
  }
  fm.out <- data.frame(matrix(0, G, length(fmM[[1]]$coef)))
  int.out <- rep(0, S)
  names(fm.out) <- names(fmM[[1]]$coef)
  tau <- data.frame(tau)
  names(tau) <- paste("grp.", 1:G, sep = "")
  EN <- -sum(unlist(tau) * log(unlist(tau)))
  d <- length(unlist(fm.out)) + length(tau) - 1
  fm.theta <- rep(0, G)
  int.out <- fmM[[1]]$sp.intercept
  for (i in 1:G) {
    fm.out[i, ] <- fmM[[i]]$coef
  }
  names(pi) <- paste("G", 1:G, sep = ".")
  t.pi <- additive_logistic(pi, TRUE)
  parms <- c(t.pi[1:(G - 1)], unlist(fm.out), int.out, rep(1,
                                                           S))
  logL.full <- logL
  logL <- logLmix_nbinom(parms, first.fit, G, S, sp, sp.name)
  options(warn = temp.warn)
  if (control$em_full_model)
    return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau,
                                                                  4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 *
                  logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out,
                theta = fm.theta, fmM = fmM, model.tau = dat.tau,
                covar = 0, aic.full = -2 * logL.full + 2 * d, bic.full = -2 *
                  logL.full + log(S) * d, pars = parms))
  return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau,
                                                                4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 *
                logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out,
              theta = fm.theta, covar = 0, aic.full = -2 * logL.full +
                2 * d, bic.full = -2 * logL.full + log(S) * d, pars = parms))
}

"fitmix_nbinom.cpp" <- function (form, datsp, sp, G=2, pars=NA, control){
  if(!is.numeric(sp)){
    sp <- as.integer(factor(sp))
  }
  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  X <- stats::model.matrix(form, data = datsp[sp==sp.name[1],])
  offset <- stats::model.frame(form, data = datsp[sp==sp.name[1],])
  offset <- stats::model.offset(offset)
  if(is.null(offset)) offset <-  rep(0,n)
  y <- datsp$obs
  if(is.na(pars[1])) {
    sp.int <- rep(0.5,S)
    sp.dispersion <- rep(1,S)
    fm <- matrix(stats::runif(ncol(X)*G,-1,1),G,ncol(X))
    pars <- c(stats::runif(G-1),unlist(fm),sp.int,sp.dispersion)
  }

  gradient <- rep(0,length(pars))
  tau <- matrix(0,S,G) ##must leave this in as defines S & G

  loglike <- try(.Call("SpeciesMix",pars,y,X,sp,tau,gradient,offset,as.integer(2),PACKAGE="ecomix"))

  calc_deriv <- function(p){
    gradient <- rep(0,length(pars))
    ll <- .Call("Calculate_Gradient",p,y,X,sp,tau,gradient,offset,as.integer(2),PACKAGE="ecomix")
    return(gradient)
  }
  r.deriv <- function(p){ logLmix_nbinom(p,list(y=y,x=stats::model.matrix(form, data = datsp)),G,S,sp,sp.name,out.tau=FALSE)}

  ##print(r.grad)
  hes <- 0
  covar <- 0
  if(control$calc.hes){
    hes <- numDeriv::jacobian(calc_deriv,pars)
    dim(hes) <- rep(length(pars),2)
    dim(hes) <- rep(length(pars),2)
    covar <- try(solve(hes))
    #rownames(hes) <- colnames(hes) <- c(paste("G.",1:(G-1),sep=""),paste("G",1:G,rep(colnames(X),each=G),sep="."))
  }
  if(!is.numeric(loglike)) loglike <- 0
  pi <- pars[1:(G-1)]
  #coef <- pars[ (G):length(pars)]
  coef <- pars[-1*(1:((G-1)))]  ## remove pi
  sp.int <- coef[(length(coef)-(2*S-1)):(length(coef)-S)]
  sp.dispersion <- coef[(length(coef)-(S-1)):length(coef)]
  fm <- coef[-1*((length(coef)-(2*S-1)):length(coef))]
  offset <- stats::model.frame(form, data = datsp)
  offset <- stats::model.offset(offset)
  if(is.null(offset)) offset <-  rep(0,nrow(datsp))

  r.logl <- logLmix_nbinom(pars,list(y=y,x=stats::model.matrix(form, data = datsp),offset=offset),G,S,sp,sp.name,out.tau=TRUE)
  print(r.logl$logl)
  pi <- additive_logistic(pi)
  names(pi) <- paste("G.",1:G,sep="")
  coef <- matrix(fm,G,ncol(X))
  rownames(coef) <- paste("G.",1:G,sep="")
  colnames(coef) <- colnames(X)

  AIC <- 2*loglike + 2*length(pars)
  BIC <- 2*loglike + log(S)*length(pars)
  ##list(logl=loglike,r.logl=r.logl$logl,pi=pi,coef=coef,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient)
  list(logl=loglike,pi=pi,coef=coef,sp.intercept=sp.int,sp.dispersion=sp.dispersion,tau=round(exp(r.logl$tau),4),aic=AIC,bic=BIC,hessian=hes,gradient=gradient,covar=covar)#,r.grad=r.grad)

}


"logLmix_nbinom" <-
  function (pars,first.fit,G,S,sp,spname,out.tau=FALSE)
  {
    tau <- matrix(0,S,G)
    ##tau,out.tau=FALSE
    if(G>1){
      fm <- pars[-1*(1:((G-1)))]  ## remove pi
      ##sp.int <- fm[(length(fm)-(S*G-1)):length(fm)]
      sp.int <- fm[(length(fm)-(2*S-1)):(length(fm)-S)]
      sp.dispersion <- fm[(length(fm)-(S-1)):length(fm)]

      ##fm <- fm[-1*((length(fm)-(S*G-1)):length(fm))]
      fm <- fm[-1*((length(fm)-(2*S-1)):length(fm))]
      pi <- pars[(1:(G-1))]
      theta <- pars[G:(G+G-1)]

      ##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
      #pi <- tau[((length(tau)-(G-2)):length(tau))]
      ##dim(sp.int) <- c(S,G)
      dim(fm) <- c(G,length(fm)/G)

      ##pi[G] <- 1-sum(pi)
      pi <- additive_logistic(pi)

    } else{
      return(0)
      fm <- tau[1:(length(pars)-1)]
      dim(fm) <- c(1,length(fm))
      pi <- 1
    }



    log.like <- 0
    for(j in 1:S){
      sel.sp <- which(sp==spname[j])
      tmp.like <- rep(0,G)
      for(i in 1:G){
        lpre <- cbind(1,first.fit$x[sel.sp,])%*%c(sp.int[j],fm[i,])+first.fit$offset[sel.sp]##,i],fm[i,])
        ##tmp.like[i] <- sum(dpois(first.fit$y[sel.sp],exp(lpre),log=TRUE))
        tmp.like[i] <- sum(stats::dnbinom(first.fit$y[sel.sp],mu=exp(lpre),size=sp.dispersion[j],log=TRUE))
      }
      # print(tmp.like)
      eps <- max(tmp.like)
      log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
      tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
    }

    if(out.tau)return(list(logl=log.like,tau=tau))
    log.like
  }

"species_mix_em_nbinom" <- function (sp.form,sp.data,covar.data,G=2,control){
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)

  fmM.out <- fitmix_nbinom(sp.form,data,sp,G,maxit,control$trace)
  if(control$em_refit>1)
    for(i in 2:control$em_refit){
      fmM <- fitmix_nbinom(sp.form,data,sp,G,maxit,control$trace)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  if(control$em_calculate_hessian){
    var <- 0
    t.pi <- additive_logistic(fmM.out$pi,TRUE)
    parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef),fmM.out$sp.intercept,rep(1,S))##unlist(fmM.out$sp.intercept))
    first.fit <- list(y=data[,1],x=stats::model.matrix(sp.form,data=data)[,-1])
    fun_est_var <- function(x){-logLmix_nbinom(x,first.fit,G,S,sp,sp.name)}
    deriv <- numDeriv::jacobian(fun_est_var,parms)
    var <- solve(numDeriv::hessian(fun_est_var, parms))
    colnames( var) <- rownames( var) <- names( parms)
    fmM.out$covar <- var
  }
  if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}

"species_mix_nbinom" <-  function (sp.form,sp.data,covar.data,G=2, pars=NA, control){
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  if(control$em_prefit | G==1){
    prefit <- species_mix_em_nbinom(sp.form,sp.data,covar.data,G,control)
    pars <- prefit$pars#c(additive_logistic(prefit$pi,TRUE)[1:(G-1)],unlist(prefit$coef))
  }
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)
  if(G==1) return(prefit)
  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
  names(data)[1] <- as.character(sp.form)[2]
  sp.form <- update.formula(sp.form,obs~-1+.)
  fmM.out <- fitmix_nbinom.cpp(sp.form,data,sp,G,pars=pars,control)
  rownames(fmM.out$tau) <- sp.name ## add names to taus
  fmM.out$se <- NA
  if(control$calculate_hessian_cpp){
    fmM.out$covar <- try(solve(fmM.out$hessian))
    if(class(fmM.out$covar)!="try-error"){
      tmp <- sqrt(diag(fmM.out$covar))
      tmp <- tmp[-1*(1:((G-1)))]
      tmp <- tmp[-1*((length(tmp)-(2*S-1)):length(tmp))]
      fmM.out$se <- matrix(tmp,G,ncol(fmM.out$coef))
      colnames(fmM.out$se) <- colnames(fmM.out$coef)
      rownames(fmM.out$se) <- rownames(fmM.out$coef)
    }
  }
  if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form
  class(fmM.out) <- c("species_mix","negbin")
  fmM.out
}

"weighted_glm_nbinom" <-  function (g,first.fit,tau,n,fmM,sp){
  dat.tau <- rep(tau[,g],each=n)
  sp.name <- unique(sp)
  X <- first.fit$x
  sp.mat <- matrix(0,dim(X)[1],length(sp.name))

  for(i in seq_along(sp.name)){
    sp.mat[sp==sp.name[i],i] <- 1
  }
  X <- cbind(sp.mat,X)
  f.mix <- glm_fit_nbinom(x=X,y=first.fit$y,offset=first.fit$offset,weights=dat.tau)
  sp.intercept <- f.mix$coef[seq_along(sp.name)]
  sp.intercept[is.na(sp.intercept)] <- 0
  return(list(coef=f.mix$coef[-1:-length(sp.name)],theta=f.mix$theta,sp.intercept=sp.intercept,fitted=f.mix$fitted))#,lpre=f.mix$linear.predictors))

}

# ##### Tweedie Functions #####
#
# "apply_glm_tweedie" <- function (i,form,datsp,tau,n){
#   dat.tau <- rep(tau[,i],each=n)
#   x <- stats::model.matrix(stats::as.formula(form),data=datsp)
#   y <- datsp$obs
#   environment(form) <- environment()
#   f.mix <- fishMod::tglm(stats::as.formula(form),data=datsp,wts=dat.tau,vcov=FALSE,residuals=FALSE,trace=0)
#   sp.int <- rep(f.mix$coef[1],dim(tau)[1])
#   return(list(coef=f.mix$coef[c(-1,-(length(f.mix$coef)-1),-(length(f.mix$coef)))],phi=f.mix$coef["phi"],p=f.mix$coef["p"],sp.intercept=sp.int))
# }
#
# "create_starting_values_tweedie" <- function (S,G,n,form,datsp,control){
#   environment(form) <- environment()
#   tau <- matrix(stats::runif(S*G),S,G)
#   tau <- (tau/rowSums(tau))
#   fmM <- list()
#   for(i in 1:G){
#     pi[i] <- sum(tau[,i])/S
#   }
#   offset <- stats::model.frame(form,datsp)
#   offset <- stats::model.offset(offset)
#   if(is.null(offset)) offset <- rep(0,length(datsp$obs))
#   fmM <- surveillance::plapply(1:G,apply_glm_tweedie,form,datsp,tau,n,.parallel=control$cores, .verbose = !control$quiet)
#   first.fit <- list(x=stats::model.matrix(stats::as.formula(form),data=datsp)[,-1],y=datsp$obs,formula=form)
#   return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
# }
#
# "create_starting_values_tweedie_kmeans" <- function(S, G, n, form, datsp, tol = 0.1, control){
#   MM <- stats::model.matrix(form, datsp)
#   offset <- stats::model.frame(form,datsp)
#   offset <- stats::model.offset(offset)
#   if(is.null(offset)) offset <- rep(0,nrow(MM))
#   first.fit <- list(x = stats::model.matrix(stats::as.formula(form), data = datsp)[,
#                                                                                    -1], y = datsp$obs, offset=offset,formula = form)
#   if(class(first.fit$x)=="numeric") {first.fit$x <- matrix(first.fit$x,length(first.fit$x),1) }#deal with model matrix returning a vector for covar == 1
#   if (tol < 0 || tol >= 1)
#     stop("Minimum Prevalence % must be between 0 and 1")
#   sp.name <- 1:S
#   sp <- rep(sp.name, each = n)
#   starting.fitem <- list(intercepts = rep(0, S), alpha = rep(0,
#                                                              S))
#   all.betas <- matrix(0, nrow = S, ncol = ncol(MM))
#   colnames(all.betas) <- colnames(MM)
#   for (j in 1:S) {
#     fit <- fishMod::tglm(form, datsp[((j - 1) * n):(j * n - 1), ],
#                          vcov = FALSE, trace = 0, p = 1.6, control)
#     starting.fitem$sp.intercepts[j] <- fit$coef[1]
#     all.betas[j, ] <- fit$coef[-length(fit$coef)]
#     starting.fitem$phi[j] = fit$coef["phi"]
#   }
#   all.betas[, 1] <- 0
#   cat("Clustering...\n")
#   fmmvnorm <- stats::kmeans(x = all.betas, centers = G, iter.max = 100,
#                             nstart = 50)
#   starting.fitem$coef <- fmmvnorm$centers
#   fmM <- list()
#   for (i in 1:G) {
#     B <- matrix(rep(fmmvnorm$centers[i, ], nrow(datsp)),
#                 nrow(datsp), ncol(fmmvnorm$centers), byrow = TRUE)
#     B[, 1] <- rep(starting.fitem$sp.intercepts, each = n)
#     fitted <- exp(rowSums(MM * B)+offset)
#     fmM[[i]] <- list(coef = fmmvnorm$centers[i, 2:ncol(fmmvnorm$centers)],
#                      phi = mean(starting.fitem$phi[fmmvnorm$cluster ==
#                                                      i]), p = 1.6, sp.intercept = starting.fitem$sp.intercepts,
#                      fitted = fitted)
#   }
#   tau <- matrix(0, S, G)
#   pi <- rep(1/G, G)
#   pi <- stats::runif(G, 0.2, 0.8)
#   pi <- pi/sum(pi)
#   est.tau <- surveillance::plapply(1:S, estimate_pi_tweedie, sp, sp.name,
#                                    datsp, fmM, pi, G, first.fit,.parallel=control$cores, .verbose = !control$quiet)
#   max.newTau <- 0.8
#   alpha <- (1 - max.newTau * G)/(max.newTau * (2 - G) - 1)
#   for (j in 1:S) {
#     newTau <- (2 * alpha * est.tau[[j]]$tau - alpha + 1)/(2 *
#                                                             alpha - alpha * G + G)
#     tau[j, ] <- newTau
#   }
#   for (i in 1:G) {
#     pi[i] <- sum(tau[, i])/S
#   }
#   return(list(pi = pi, fmM = fmM, tau = tau, first.fit = first.fit))
# }
#
# "estimate_pi_tweedie" <- function (j, sp, spname, datsp, fmM, pi, G, first.fit){
#   tmp.like <- rep(0, G)
#   tau <- rep(0, G)
#   sel.sp <- which(sp == spname[j])
#   for (i in 1:G) {
#     lpre <- dTweedie(first.fit$y[sel.sp], mu = fmM[[i]]$fitted[sel.sp],
#                      phi = fmM[[i]]$phi, p = fmM[[i]]$p, LOG = TRUE)
#     tmp.like[i] <- sum(lpre)
#   }
#   eps <- max(tmp.like)
#   sum.like <- (log(sum(pi * exp((tmp.like) - (eps)))) + (eps))
#   for (i in 1:G) {
#     tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like)
#   }
#   return(list(tau = tau, sum.like = sum.like))
# }
#
#
#
# "fitmix_tweedie" <-  function (form, datsp, sp, G = 2, control){
#   temp.warn <- getOption("warn")
#   options(warn = -1)
#   sp.name <- unique(sp)
#   S <- length(unique(sp))
#   n <- length(which(sp == sp.name[1]))
#   cat("Fitting Group", G, "\n")
#   if (trace)
#     cat("Iteration | LogL \n")
#   dat.tau <- 0
#   pi <- rep(0, G)
#   ite <- 1
#   logL <- -99999999
#   old.logL <- -88888888
#   if (ite != 1)
#     t1 <- create_starting_values_tweedie(S, G, n, form, datsp)
#   t1 <- create_starting_values_tweedie_kmeans(S, G, n, form,
#                                               datsp)
#   pi <- t1$pi
#   fmM <- t1$fmM
#   tau <- t1$tau
#   dat.tau <- t1$dat.tau
#   first.fit <- t1$first.fit
#   while (abs(logL - old.logL) > 1e-04 & ite <= maxit) {
#     old.logL <- logL
#     for (i in 1:G) {
#       pi[i] <- sum(tau[, i])/S
#     }
#     if (any(pi == 0)) {
#       cat("pi has gone to zero - restarting fitting \n")
#       t1 <- create_starting_values_tweedie(S, G, n, form,
#                                            datsp)
#       pi <- t1$pi
#       fmM <- t1$fmM
#       tau <- t1$tau
#       dat.tau <- t1$dat.tau
#       first.fit <- t1$first.fit
#       ite <- 1
#     }
#     fmM <- surveillance::plapply(1:G, weighted_glm_tweedie, first.fit, tau, n, fmM, sp, .parallel=control$cores, .verbose = !control$quiet)
#     for (j in 1:S) {
#       tmp <- rep(0, G)
#       for (g in 1:G) tmp[g] <- fmM[[g]]$sp.intercept[j]
#       tmp <- sum(tmp * tau[j, ])
#       for (g in 1:G) fmM[[g]]$sp.intercept[j] <- tmp
#     }
#     logL <- 0
#     tmp.like <- matrix(0, S, G)
#     est.tau <- surveillance::plapply(1:S, estimate_pi_tweedie, sp, sp.name,
#                                      datsp, fmM, pi, G, first.fit, .parallel=control$cores, .verbose = !control$quiet)
#     for (j in 1:S) {
#       if (is.atomic(est.tau[[j]])) {
#         print(est.tau[[j]])
#       }
#       else {
#         tau[j, ] <- est.tau[[j]]$tau
#         logL <- logL + est.tau[[j]]$sum.like
#       }
#     }
#     if (trace)
#       cat(ite, "     | ", logL, "\n")
#     ite <- ite + 1
#   }
#   fm.out <- data.frame(matrix(0, G, length(fmM[[1]]$coef)))
#   int.out <- fmM[[1]]$sp.intercept
#   names(fm.out) <- names(fmM[[1]]$coef)
#   tau <- data.frame(tau)
#   names(tau) <- paste("grp.", 1:G, sep = "")
#   EN <- -sum(unlist(tau) * log(unlist(tau)))
#   d <- length(unlist(fm.out)) + length(tau) - 1
#   fm.phi <- fm.p <- rep(0, G)
#   for (i in 1:G) {
#     fm.out[i, ] <- fmM[[i]]$coef
#     fm.phi[i] <- fmM[[i]]$phi
#     fm.p[i] <- fmM[[i]]$p
#   }
#   names(pi) <- paste("G", 1:G, sep = ".")
#   t.pi <- additive_logistic(pi, TRUE)
#   parms <- c(t.pi[1:(G - 1)], unlist(fm.out), int.out, rep(2,
#                                                            S))
#   names(parms) <- c(rep("pi", G - 1), rep("coef", length(unlist(fm.out))),
#                     rep("int", length(int.out)), rep("phi", S))
#   logL.full <- logL
#   logL <- logLmix_tweedie(parms, first.fit, G, S, sp, sp.name)
#   print(logL)
#   options(warn = temp.warn)
#   if (control$em_full_model)
#     return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau,
#                                                                   4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 *
#                   logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out,
#                 phi = fm.phi, p = fm.p, fmM = fmM, model.tau = dat.tau,
#                 covar = 0, aic.full = -2 * logL.full + 2 * d, bic.full = -2 *
#                   logL.full + log(S) * d, pars = parms))
#   return(list(logl = logL, aic = -2 * logL + 2 * d, tau = round(tau,
#                                                                 4), pi = pi, bic = -2 * logL + log(S) * d, ICL = -2 *
#                 logL + log(S) * d + 2 * EN, coef = fm.out, sp.intercept = int.out,
#               phi = fm.phi, p = fm.p, covar = 0, aic.full = -2 * logL.full +
#                 2 * d, bic.full = -2 * logL.full + log(S) * d, pars = parms))
# }
#
#
# "fitmix_tweedie.cpp" <-
#   function (form, datsp, sp, G = 2, pars = NA, trace = TRUE, calc.hes = FALSE)
#   {
#     if (!is.numeric(sp)) {
#       sp <- as.integer(factor(sp))
#     }
#     sp.name <- unique(sp)
#     S <- length(unique(sp))
#     n <- length(which(sp == sp.name[1]))
#     X <- stats::model.matrix(form, data = datsp[sp == sp.name[1], ])
#     offset <- stats::model.frame(form, data = datsp[sp==sp.name[1],])
#     offset <- stats::model.offset(offset)
#     if(is.null(offset)) offset <-  rep(0,n)
#
#     y <- datsp$obs
#     if (is.na(pars[1])) {
#       sp.int <- rep(0.5, S)
#       sp.phi <- rep(1, S)
#       sp.p <- rep(1.6, S)
#       fm <- matrix(rep(0.5, ncol(X) * G), G, ncol(X))
#       pars <- c(rep(0.5, G - 1), unlist(fm), sp.int, sp.phi)
#     }
#     pars.og <- pars
#     pars.og[1] <- pars.og[1] * 0.99
#     gradient <- rep(0, length(pars))
#     tau <- matrix(0, S, G)
#     loglike <- try(.Call("SpeciesMix", pars, y, X, sp, tau, gradient,
#                          offset, as.integer(3),PACKAGE="ecomix"))
#     calc_deriv <- function(p) {
#       gradient <- rep(0, length(pars))
#       ll <- .Call("Calculate_Gradient", p, y, X, sp, tau, gradient,
#                   offset, as.integer(3),PACKAGE="ecomix")
#       return(gradient)
#     }
#     hes <- 0
#     if (calc.hes) {
#       hes <- numDeriv::jacobian(calc_deriv, pars, method = "simple")
#       dim(hes) <- rep(length(pars), 2)
#       dim(hes) <- rep(length(pars), 2)
#     }
#     if (!is.numeric(loglike))
#       loglike <- 0
#     pi <- pars[1:(G - 1)]
#     coef <- pars[-1 * (1:((G - 1)))]
#     sp.int <- coef[(length(coef) - (2 * S - 1)):(length(coef) -
#                                                    S)]
#     sp.dispersion <- coef[(length(coef) - (S - 1)):length(coef)]
#     sp.phi <- sp.dispersion[1:S]
#     sp.p <- rep(1.6, S)
#     fm <- coef[-1 * ((length(coef) - (2 * S - 1)):length(coef))]
#     offset <- stats::model.frame(form, data = datsp)
#     offset <- stats::model.offset(offset)
#     if(is.null(offset)) offset <- rep(0,nrow(datsp))
#
#     names(gradient) <- names(pars) <- c(rep("pi", G - 1), rep("coef",
#                                                               length(unlist(fm))), rep("int", length(sp.int)), rep("phi",
#                                                                                                                    S))
#     r.deriv <- function(p) {
#       logLmix_tweedie(p, list(y = y, x = stats::model.matrix(form,
#                                                              data = datsp)), G, S, sp, sp.name, out.tau = FALSE)
#     }
#     r.logl <- logLmix_tweedie(pars, list(y = y, x = stats::model.matrix(form,
#                                                                         data = datsp),offset=offset), G, S, sp, sp.name, out.tau = TRUE)
#     print(r.logl$logl)
#     pi <- additive_logistic(pi)
#     names(pi) <- paste("G.", 1:G, sep = "")
#     coef <- matrix(fm, G, ncol(X))
#     rownames(coef) <- paste("G.", 1:G, sep = "")
#     colnames(coef) <- colnames(X)
#     AIC <- 2 * loglike + 2 * length(pars)
#     BIC <- 2 * loglike + log(S) * length(pars)
#     list(logl = loglike, pi = pi, coef = coef, sp.intercept = sp.int,
#          phi = sp.phi, p = sp.p, tau = round(exp(r.logl$tau),
#                                              4), aic = AIC, bic = BIC, hessian = hes, gradient = gradient)
#   }
# "ldTweedie.lp" <- function ( parms, y, X.p, offsetty, phi, p, wts=rep( 1, length( y))){
#   mu <- exp( X.p %*% parms[seq_len(ncol(X.p))] + offsetty)
#
#   if( is.null( phi) & is.null( p)){
#     phi <- parms[ncol( X.p) + 1]
#     p <- parms[ncol( X.p) + 2]
#   }
#   if( is.null( phi) & !is.null( p))
#     phi <- parms[ncol( X.p)+1]
#   if( !is.null( phi) & is.null( p))
#     p <- parms[ncol( X.p)+1]
#
#   lambda <- ( mu^( 2-p)) / ( phi*(2-p))
#   alpha <- ( 2-p) / ( p-1)
#   tau <- phi*(p-1)*mu^(p-1)
#   mu.Z <- alpha * tau
#
#   return( -sum( wts * dPoisGam( y, lambda=lambda, mu.Z=mu.Z, alpha=alpha, LOG=TRUE)))
# }
#
#
# "ldTweedie.lp.deriv" <-
#   function ( parms, y, X.p, offsetty, phi, p, wts=rep( 1, length( y)))
#   {
#     mu <- exp( X.p %*% parms[seq_len(ncol( X.p))] + offsetty)
#
#     p.flag <- phi.flag <- FALSE
#     if( is.null( phi) & is.null( p)){
#       p.flag <- phi.flag <- TRUE
#       phi <- parms[ncol( X.p) + 1]
#       p <- parms[ncol( X.p) + 2]
#     }
#     if( is.null( phi) & !is.null( p)){
#       phi <- parms[ncol( X.p)+1]
#       phi.flag <- TRUE
#     }
#     if( !is.null( phi) & is.null( p)){
#       p <- parms[ncol( X.p)+1]
#       p.flag <- TRUE
#     }
#
#     lambda <- ( mu^( 2-p)) / ( phi*(2-p))
#     alpha <- ( 2-p) / ( p-1)
#     tau <- phi*(p-1)*mu^(p-1)
#     mu.Z <- alpha * tau
#
#     dTweedparms <- -wts * dPoisGamDerivs( y, lambda=lambda, mu.Z=mu.Z, alpha=alpha)
#
#     DTweedparmsDmu <- matrix( c( ( mu^(1-p)) / phi, alpha*phi*( ( p-1)^2)*( mu^(p-2)), rep( 0, length( mu))), nrow=3, byrow=TRUE)
#     tmp <- rowSums( dTweedparms * t( DTweedparmsDmu))
#     tmp <- tmp * mu
#     tmp <- apply( X.p, 2, function( x) x*tmp)
#
#     derivs <- colSums( tmp)
#
#     if( phi.flag){
#       DTweedparmsDphi <- matrix( c( -( ( mu^(2-p)) / ( ( phi^2)*(2-p))), alpha*( p-1)*( mu^( p-1)), rep( 0, length( mu))), nrow=3, byrow=TRUE)
#       tmpPhi <- rowSums( dTweedparms * t( DTweedparmsDphi))#vectorised way of doing odd calculation
#       derivs <- c( derivs, sum( tmpPhi))
#       names( derivs)[length( derivs)] <- "phi"
#     }
#     if( p.flag){
#       dalphadp <- -( 1+alpha) / ( p-1)
#       DTweedparmsDp <- matrix( c( lambda*( 1/(2-p) - log( mu)), mu.Z*( dalphadp/alpha + 1/( p-1) + log( mu)), rep( dalphadp, length( y))), nrow=3, byrow=TRUE)
#       tmpP <- rowSums( dTweedparms * t( DTweedparmsDp))
#       derivs <- c( derivs, sum( tmpP))
#       names( derivs)[length( derivs)] <- "p"
#     }
#
#     return( derivs)
#   }
#
# "pTweedie" <-  function ( quant, mu, phi, p){
#   lambda <- ( mu^( 2-p)) / ( phi*(2-p))
#   alpha <- ( 2-p) / ( p-1)
#   tau <- phi*(p-1)*mu^(p-1)
#   mu.Z <- alpha * tau
#   ps <- 0
#   return( ps)
# }
#
# "rTweedie" <- function ( n, mu, phi, p) {
#   lambda <- ( mu^( 2-p)) / ( phi*(2-p))
#   alpha <- ( 2-p) / ( p-1)
#   tau <- phi*(p-1)*mu^(p-1)
#   mu.Z <- alpha * tau
#   rans <- rPoisGam( n, lambda, mu.Z, alpha)
#   return( rans)
# }
#
# "species_mix_em_tweedie" <- function (sp.form,sp.data,covar.data,G=2,control){
#   S <- dim(sp.data)[2]
#   if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
#   n <- dim(sp.data)[1]
#   sp <- rep(sp.name,each=n)
#
#   var.names <- colnames(covar.data)
#   covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
#   names(covar.data) <- var.names
#   sp.data <- as.data.frame(sp.data)
#   data <- data.frame(obs=unlist(sp.data),covar.data)
#
#   fmM.out <- fitmix_tweedie(sp.form,data,sp,G,maxit,control$trace)
#   if(control$em_refit>1)
#     for(i in 2:control$em_refit){
#       fmM <- fitmix_tweedie(sp.form,data,sp,G,maxit,control$trace)
#       if(fmM$logl>fmM.out$logl) fmM.out <- fmM
#     }
#   # control$em_calculate_hessian <- FALSE
#   if(control$em_calculate_hessian){
#     var <- 0
#     t.pi <- additive_logistic(fmM.out$pi,TRUE)
#     parms <- c(t.pi[1:(G-1)],unlist(fmM.out$coef),unlist(fmM.out$sp.intercept),fmM.out$theta)
#     first.fit <- list(y=data[,1],x=stats::model.matrix(sp.form,data=data)[,-1])
#     fun_est_var <- function(x){-logLmix_tweedie(x,first.fit,G,S,sp,sp.name)}
#     var <- solve(numDeriv::hessian(fun_est_var, parms))
#     colnames( var) <- rownames( var) <- names( parms)
#     fmM.out$covar <- var
#   }
#   if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
#   fmM.out$formula <- sp.form
#
#   fmM.out
# }
#
# "species_mix_tweedie" <- function (sp.form,sp.data,covar.data,G=2, pars=NA, em_prefit=TRUE,control){
#   t.covar.data <- covar.data
#   t.sp.data <- sp.data
#   sp.form <- update.formula(sp.form,obs~1+.)
#   pars <- NA
#   if(control$em_prefit | G==1){
#     prefit <- species_mix_em_tweedie(sp.form,sp.data,covar.data,G,control)
#     pars <- prefit$pars##pars <- c(additive_logistic(prefit$pi,TRUE)[1:(G-1)],unlist(prefit$coef))
#   }
#   S <- dim(sp.data)[2]
#   if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
#   n <- dim(sp.data)[1]
#   sp <- rep(sp.name,each=n)
#   if(G==1) return(prefit)
#   var.names <- colnames(covar.data)
#   covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
#   names(covar.data) <- var.names
#   sp.data <- as.data.frame(sp.data)
#   data <- data.frame(obs=as.numeric(unlist(sp.data)),covar.data)
#   names(data)[1] <- as.character(sp.form)[2]
#   sp.form <- update.formula(sp.form,obs~-1+.)
#   fmM.out <- fitmix_tweedie.cpp(sp.form,data,sp,G,pars=pars,calc.hes=control$calculate_hessian)
#   rownames(fmM.out$tau) <- sp.name ## add names to taus
#   fmM.out$se <- NA
#   if(control$calculate_hessian){
#     fmM.out$covar <- try(solve(fmM.out$hessian))
#     if(class(fmM.out$covar)!="try-error"){
#       colnames(fmM.out$covar) <- rownames(fmM.out$covar) <- names(fmM.out$gradient)
#       tmp <- sqrt(diag(fmM.out$covar))
#       tmp <- tmp[-1*(1:((G-1)))]
#       tmp <- tmp[-1*((length(tmp)-(2*S-1)):length(tmp))]
#       fmM.out$se <- matrix(tmp,G,ncol(fmM.out$coef))
#       colnames(fmM.out$se) <- colnames(fmM.out$coef)
#       rownames(fmM.out$se) <- rownames(fmM.out$coef)
#     }
#   }
#   if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
#   fmM.out$formula <- sp.form
#   class(fmM.out) <- c("species_mix","tweedie")
#   fmM.out
# }
#
#
#
# "weighted_glm_tweedie" <-  function (g, first.fit, tau, n, fmM, sp, control){
#   dat.tau <- rep(tau[, g], each = n)
#   sp.int.offset <- rep(fmM[[g]]$sp.intercept, each = n)+first.fit$offset
#   sp.name <- unique(sp)
#   X <- first.fit$x
#   f.mix <- fishMod::tglm.fit(x = X, y = first.fit$y, wts = dat.tau,
#                              offset = sp.int.offset, control=control, inits = fmM[[g]]$coef, phi = fmM[[g]]$phi,
#                              p = 1.6)
#   return(list(coef = f.mix$coef, phi = fmM[[g]]$phi, p = 1.6,
#               sp.intercept = fmM[[g]]$sp.intercept, fitted = f.mix$fitted))
# }

##### Gaussian functions #####

"apply_glm_gaussian" <-  function (i,form,datsp,tau,n){
  dat.tau <- rep(tau[,i],each=n)
  x <- stats::model.matrix(stats::as.formula(form),data=datsp)
  y <- datsp$obs
  environment(form) <- environment()
  f.mix <- stats::glm(stats::as.formula(form),data=datsp,weights=dat.tau,family="gaussian")
  sp.int <- rep(f.mix$coef[1],dim(tau)[1])
  return(list(coef=f.mix$coef[-1],theta=sqrt(f.mix$deviance/f.mix$df.residual),sp.intercept=sp.int))
}

"create_starting_values_gaussian" <- function (S,G,n,form,datsp,control){
  environment(form) <- environment()
  tau <- matrix(stats::runif(S*G),S,G)
  tau <- (tau/rowSums(tau))
  fmM <- list()
  for(i in 1:G){
    pi[i] <- sum(tau[,i])/S
  }
  fmM <- surveillance::plapply(1:G,apply_glm_gaussian,form,datsp,tau,n,.parallel=control$cores, .verbose = !control$quiet)
  first.fit <- list(x=stats::model.matrix(stats::as.formula(form),data=datsp)[,-1],y=datsp$obs,formula=form)
  return(list(pi=pi,fmM=fmM,tau=tau,first.fit=first.fit))
}

"estimate_pi_gaussian" <-  function (j,sp,spname,datsp,fmM,pi,G,first.fit){

  tmp.like <- rep(0,G)
  tau <- rep(0,G)
  sel.sp <- which(sp==spname[j])

  for(i in 1:G) {
    lpre <- stats::dnorm(first.fit$y[sel.sp],mean=fmM[[i]]$fitted[sel.sp],sd=1,log=TRUE)
    tmp.like[i] <- sum(lpre)
  }

  eps <- max(tmp.like)
  sum.like <- (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
  for(i in 1:G) {
    tau[i] <- exp((log(pi[i]) + tmp.like[i]) - sum.like )
  }
  return(list(tau=tau,sum.like=sum.like))
}

"fitmix_gaussian" <- function (form,datsp,sp,G=2,maxit=500,control){
  ## fitting abundance data with a negative binomial distribution
  ## dat2 has colums obs,sp
  ##
  temp.warn <- getOption( "warn")
  options( warn=-1)

  sp.name <- unique(sp)
  S <- length(unique(sp))
  n <- length(which(sp==sp.name[1]))

  cat("Fitting Group",G,"\n")
  if(control$trace==1) cat("Iteration | LogL \n")

  ##dat.tau <- data.frame(matrix(0,dim(datsp)[1],G))
  dat.tau <- 0
  ##dat <- data.frame(datsp,dat.tau)
  pi <- rep(0,G)
  ite <- 1
  logL <- -99999999
  old.logL <- -88888888

  ## set up initial GLM
  t1 <- create_starting_values_gaussian(S,G,n,form,datsp)
  pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit


  while(abs(logL-old.logL) > 0.0001 & ite<=maxit){
    old.logL <- logL

    for(i in 1:G){
      ##dat.tau[,i] <- rep(tau[,i],each=n)
      pi[i] <- sum(tau[,i])/S
    }

    if(any(pi==0)) { ## occasionally with complicated models the random starting values result in a pi[i]==0; so restart with new random starts
      cat("pi has gone to zero - restarting fitting \n")
      t1 <- create_starting_values_gaussian(S,G,n,form,datsp)
      pi <- t1$pi;fmM <- t1$fmM;tau <- t1$tau;dat.tau <- t1$dat.tau;first.fit <- t1$first.fit
      ite <- 1
    }

    fmM <- surveillance::plapply(1:G,weighted_glm_gaussian,first.fit,tau,n,fmM,sp,.parallel=cores, .verbose = !control$quiet)


    logL <- 0
    tmp.like <- matrix(0,S,G)

    est.tau <- surveillance::plapply(1:S,estimate_pi_gaussian,sp,sp.name,datsp,fmM,pi,G,first.fit, .parallel=control$cores, .verbose = !control$quiet)

    for(j in 1:S){
      if(is.atomic(est.tau[[j]])){ print (est.tau[[j]])} else
      {
        tau[j,] <- est.tau[[j]]$tau
        logL <- logL+est.tau[[j]]$sum.like
      }
    }

    if(control$trace==1) cat(ite,"     | ",logL,"\n")
    ite <- ite+1
  }
  fm.out <- data.frame(matrix(0,G,length(fmM[[1]]$coef)))
  int.out <- data.frame(matrix(0,S,G))
  names(int.out) <- paste("grp.",1:G,sep="")
  names(fm.out) <- names(fmM[[1]]$coef)
  tau <- data.frame(tau)
  names(tau) <- paste("grp.",1:G,sep="")
  EN <- -sum(unlist(tau)*log(unlist(tau)))
  d <- length(unlist(fm.out)) + length(tau)-1
  fm.theta <- rep(0,G)
  for(i in 1:G) {
    fm.out[i,] <- fmM[[i]]$coef
    int.out[,i] <- fmM[[i]]$sp.intercept
    ##dat.tau[,i] <- rep(tau[,i],each=n)
    fm.theta[i] <- fmM[[i]]$theta
  }

  names(pi) <- paste("G",1:G,sep=".")
  t.pi <- additive_logistic(pi,TRUE)
  parms <- c(t.pi[1:(G-1)],fm.theta,unlist(fm.out),unlist(int.out))
  logL.full <- logL
  logL <- logLmix_gaussian(parms,first.fit,G,S,sp,sp.name)

  options(warn=temp.warn)
  if(control$em_full_model)  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,sp.intercept=int.out,theta=fm.theta,fmM=fmM,model.tau=dat.tau,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))

  return(list(logl=logL,aic = -2*logL + 2*d,tau=round(tau,4),pi=pi,bic=-2*logL + log(S)*d,ICL= -2*logL + log(S)*d +2*EN,coef=fm.out,sp.intercept=int.out,theta=fm.theta,covar=0,aic.full= -2*logL.full + 2*d,bic.full= -2*logL.full + log(S)*d))

}

"logLmix_gaussian" <-
  function (pars,first.fit,G,S,sp,spname,out.tau=FALSE)
  {
    tau <- matrix(0,S,G)
    ##tau,out.tau=FALSE
    if(G>1){
      fm <- pars[-1*(1:((G-1)+G))]  ## remove pi
      sp.int <- fm[(length(fm)-(S*G-1)):length(fm)]
      fm <- fm[-1*((length(fm)-(S*G-1)):length(fm))]
      pi <- pars[(1:(G-1))]
      theta <- pars[G:(G+G-1)]

      ##    fm <- tau[-1*((length(tau)-(G-2)):length(tau))]
      #pi <- tau[((length(tau)-(G-2)):length(tau))]
      dim(sp.int) <- c(S,G)
      dim(fm) <- c(G,length(fm)/G)
      ##pi[G] <- 1-sum(pi)
      pi <- additive_logistic(pi)

    } else{
      return(0)
      fm <- tau[1:(length(pars)-1)]
      dim(fm) <- c(1,length(fm))
      pi <- 1
    }



    log.like <- 0
    for(j in 1:S){
      sel.sp <- which(sp==spname[j])
      tmp.like <- rep(0,G)
      for(i in 1:G){
        lpre <- cbind(1,first.fit$x[sel.sp,])%*%c(sp.int[j,i],fm[i,])
        ##tmp.like[i] <- sum(dpois(first.fit$y[sel.sp],exp(lpre),log=TRUE))
        tmp.like[i] <- sum(stats::dnorm(first.fit$y[sel.sp],mean=lpre,sd=theta[i],log=TRUE))
        eps <- max(tmp.like)
        log.like <- log.like +  (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
        tau[j,] <- log(pi) + tmp.like - (log(sum(pi*exp((tmp.like)-(eps))))+(eps))
      }
    }
    if(out.tau)return(list(logl=log.like,tau=tau))
    log.like
  }

"species_mix_em_gaussian" <- function (sp.form,sp.data,covar.data,G=2,control){
  S <- dim(sp.data)[2]
  if(is.null(colnames(sp.data))){sp.name <- 1:S} else {sp.name <- colnames(sp.data)}
  n <- dim(sp.data)[1]
  sp <- rep(sp.name,each=n)

  var.names <- colnames(covar.data)
  covar.data <- data.frame(kronecker( rep( 1, S), as.matrix(covar.data)))
  names(covar.data) <- var.names
  sp.data <- as.data.frame(sp.data)
  data <- data.frame(obs=unlist(sp.data),covar.data)

  fmM.out <- fitmix_gaussian(sp.form,data,sp,G,control$maxit,control$trace)
  if(control$em_refit>1)
    for(i in 2:control$em_refit){
      fmM <- fitmix_gaussian(sp.form,data,sp,G,control$maxit,control$trace)
      if(fmM$logl>fmM.out$logl) fmM.out <- fmM
    }
  if(control$em_calculate_hessian){
    var <- 0
    t.pi <- additive_logistic(fmM.out$pi,TRUE)
    parms <- c(t.pi[1:(G-1)],fmM.out$theta,unlist(fmM.out$coef),unlist(fmM.out$sp.intercept))
    first.fit <- list(y=data[,1],x=stats::model.matrix(sp.form,data=data)[,-1])
    fun_est_var <- function(x){-logLmix_gaussian(x,first.fit,G,S,sp,sp.name)}
    var <- solve(numDeriv::hessian(fun_est_var, parms))
    colnames( var) <- rownames( var) <- names( parms)
    fmM.out$covar <- var
  }
  # residuals <- FALSE
  if(control$residuals) fmM.out$residuals <- mix.residuals(fmM.out,sp.form,data,sp)
  fmM.out$formula <- sp.form

  fmM.out
}

"species_mix_gaussian" <- function (sp.form,sp.data,covar.data,G=2, pars=NA, control){
  t.covar.data <- covar.data
  t.sp.data <- sp.data
  sp.form <- update.formula(sp.form,obs~1+.)
  control$em_steps=100
  prefit <- species_mix_em_gaussian(sp.form,sp.data,covar.data,G,control)
  return(prefit)
}
"weighted_glm_gaussian" <-  function (g,first.fit,tau,n,fmM,sp){
  dat.tau <- rep(tau[,g],each=n)
  sp.name <- unique(sp)
  X <- first.fit$x
  sp.mat <- matrix(0,dim(X)[1],length(sp.name))

  for(i in seq_along(sp.name)){
    sp.mat[sp==sp.name[i],i] <- 1
  }
  X <- cbind(sp.mat,X)
  f.mix <- stats::glm.fit(x=X,y=first.fit$y,weights=dat.tau,family=gaussian())
  sp.intercept <- f.mix$coef[seq_along(sp.name)]
  sp.intercept[is.na(sp.intercept)] <- 0
  return(list(coef=f.mix$coef[-1:-length(sp.name)],theta=sqrt(f.mix$deviance/f.mix$df.residual),sp.intercept=sp.intercept,fitted=f.mix$fitted))#,lpre=f.mix$linear.predictors))

}
