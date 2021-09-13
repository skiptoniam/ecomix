##### Main species mix functions to export #####
#' @title This is how you fit a species archetype model (SAM) in ecomix.
#' @rdname species_mix
#' @name species_mix
#' @description Fits a mixture-of-regressions (SAM) to identify species archetypes.
#' @details species_mix is used to fit mixtures of glms to multivariate
#' species data. The function uses BFGS to optimize the mixture likelihood.
#' There is the option to use ECM algorithm to get appropriate starting
#' parameters. `species_mix` acts as a wrapper for species_mix.fit
#' that allows for easier data input. The data frames are merged into
#' the appropriate format for the use in species_mix.fit.
#' Minima is found using vmmin (BFGS). Currently 'bernoulli', 'binomial',
#' 'poisson', 'negative.binomial' and 'gaussian' distributions can be fitted
#' using the species_mix function. For Point Process extensions of SAMs please
#' refer to the \link{ecomix.ppm} package.
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
#'  you might use this an alternative to an offset, where there might be a set
#'  covariates which means artifact of how the data was collected or some other
#'  process like a seasonality.
#' @param data a matrix or data.frame which contains the 'species_data'
#' matrix, a const and the covariates in the structure of spp1, spp2, spp3,
#' const, temperature, rainfall. dims of matrix should be
#' nsites*(nspecies+const+covariates).
#' @param nArchetypes The number of archetypes (mixing components/groups) to
#' estimate from the data. This need to be explicitly declared. By default it is
#' set to 3.
#' @param family The family of statistical family to use within
#' the ecomix models. family can be a function from \link[stats]{family}, such as
#' binomial(). OR it can be a choice between "bernoulli", "binomial", "poisson",
#' "negative.binomial", "tweedie" and "gaussian" families. An error will be
#' thrown if the family is not recognized. The family you use is generally
#' specific to the type of data you have. One can also specify a link function
#' such as "cloglog" via binomial(link="cloglog"). If you wish to fit a point
#' process version of a SAM please refer to the \link{ecomix.ppm} package.
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
#' @param power The power parameter for 'tweedie' family. Default is 1.6, and
#' this is assigned to all species
#' @param control a list of control parameters for optimization and calculation.
#' See details.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficients. These are family specific, but at a
#' minimum you will need pi (additive_logistic transformed), alpha
#' (intercepts) and beta (mixing coefs).
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
#' @details species_mix is used to fit mixtures of glms to multivariate
#' species data. The function uses BFGS to optimize the mixture likelihood.
#' There is the option to use ECM algorithm to get appropriate starting
#' parameters. `species_mix` acts as a wrapper for species_mix.fit
#' that allows for easier data input. The data frames are merged into
#' the appropriate format for the use in species_mix.fit.
#' Minima is found using vmmin (BFGS).
#'
#' Control arguments for species_mix are as follows:
#' \describe{
#'  \item{maxit}{The maximum number of iterations. Default is 1000.}
#'  \item{quiet}{Should any reporting be performed? Default is FALSE, for reporting.}
#'  \item{trace}{Non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information. If quiet=TRUE, then trace will be set to zero.}
#'  \item{init.glmnet}{Use glmnet to get initial starting values, default is TRUE, if FALSE, glm2 will be used.}
#'  \item{init.method}{Which initialsation approach to use for covariates, default is "random2" which is uses k-means on the groups and adds some random noise. Alternatives are 'kmeans' and 'kmed'}
#'  \item{init.sd}{The amount of gaussian variation to add to starting values, default is NA, which adds variantion in based on the standard deviation of the parameters. Otherwise should be a double like 0.3}
#'  \item{minimum.presencess}{The number of presences required for each species to be included when generaling the starting values, the default is to remove all species <= 30 occurrences during inititalisation, these species will still be included in the model fit.}
#'  \item{ecm.prefit}{Use ECM to generate starting values? If FALSE, the model will use starting values without the ECM.}
#'  \item{ecm.refits}{The number of ECM refits to run, the default is 2, users could increase this to potentially get better starting values.}
#'  \item{ecm.steps}{The numer of ECM steps/iterations to run per ECM fit, the default is 5.}
#'  \item{ecm.reltol}{The relative tolerance for the ECM, default is quiet large 1e-2, makin this smaller will make the convergance tolerance stickter, you'll have to also adjust the "ecm.steps", to be larger if you make this stickter.}
#'  \item{print.cpp.inits}{Print the intialisation values being passed to c++. Handy for debugging especially if 'vmmin' error call}
#'  \item{nreport}{The frequency of reports for optimisation. Default is 10 a report for 10th iteration.}
#'  \item{abstol}{The abstol for vmmin, the default is -Inf. See \link[stats]{optim} for further details}
#'  \item{reltol}{Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.}
#'  \item{optimise.cpp}{Should optimisation for estimation occur? If TRUE (default) optimisation will occur. If FALSE no optimisation is performed.}
#'  \item{getscores.cpp}{Should we return the scores (derivates) when doing optimsation? If FALSE (default) scores will not be returned. If TRUE scores will be returned.}
#'  \item{loglOnly.cpp}{Should the log-likelihood be caulcated? If TRUE (default) then log-likelihood is calculated and returned. If FALSE then the log-likelihood is not calculated for return.}
#'  \item{derivOnly.cpp}{Should the scores be evaluated at the (final) parameter values. If TRUE (default) then they are calculated. If FALSE then they are not calculated.}
#'  \item{doPenalties}{A boolean 0 or 1. Default is 0 (false) and no penalities will be applied within the optimsation.}
#'  \item{penalty.pi}{A numeric scalar. This is the penalty for the mixing coefs. The penality is an from a Dirichlet distribution.}
#'  \item{penalty.alpha}{A numeric scalar. This is the penalty for the alpha parameters in the species model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
#'  \item{penalty.beta}{A numeric scalar. This is the penalty for the beta parameters in the  group model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
#'  \item{penalty.gamma}{A numeric scalar. This is the penalty for the gamma parameters in the species model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
#'  \item{penalty.delta}{A numeric scalar. This is the penalty for the delta parameters in the bias model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
#'  }
#'
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
                          data,
                          nArchetypes = 3,
                          family="bernoulli", offset=NULL,
                          weights=NULL, bb_weights=NULL, size = NULL, power=1.6,
                          control=list(), inits=NULL, titbits = TRUE){

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

  # setup the model.frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.exclude"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))
  dat <- clean_data_sam(mf, archetype_formula, species_formula, all_formula)#, family)

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
    message( "There are ", nArchetypes, " archetypes")

  # get archetype model matrix
  Xdat <- get_X_sam(mf.X = dat$mf.X)
  X <- Xdat$X
  xterms <- Xdat$mt.x

  # what is the W matrix (species covariates)
  Wdat <- get_W_sam(mf.W = dat$mf.W)
  W <- Wdat$W
  wterms <- Wdat$mt.w

  # what is the U matrix (all covariates)
  Udat <- get_U_sam(mf.U = dat$mf.U)
  U <- Xdat$U
  uterms <- Xdat$mt.u

  # collect terms
  tt <- list(xterms=xterms,wterms=wterms,uterms=uterms)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names)#,family)
  spp_weights <- check_spp_weights(bb_weights,S)

  # get size for binomial
  size <- check_size_binomial(size,nrow(dat$mf.X))

  # get family info
  disty.listy <- get_family_sam(family,size)
  disty <- disty.listy$disty
  linky <- disty.listy$link
  family <- disty.listy$family

  # check powers
  powers <- get_power_sam(disty,power,S)

  # summarizing data to console
  if(!control$quiet) print_input_sam(y, X, W, U, S, archetype_formula,
                                     species_formula, all_formula,
                                     family, linky,
                                     quiet=control$quiet)

  # fit species mix.
  G <- nArchetypes
  tmp <- species_mix.fit(y=y, X=X, W=W, U=U, G=G, S=S,
                         spp_weights=spp_weights,
                         site_spp_weights=site_spp_weights,
                         offset=offset, disty=disty, linky=linky, y_is_na=y_is_na,
                         size=size, powers=powers, control=control, inits=inits)

  tmp$family <- family
  tmp$disty <- disty
  tmp$link <- linky

  tmp$S <- S;
  tmp$G <- G;
  tmp$npx <- ncol(X);
  tmp$npw <- ifelse(ncol(W)>1,ncol(W),0);
  tmp$npu <- ifelse(!is.null(U),ncol(U),0);

  if(nArchetypes==1){
    tmp$pi <- tmp$pi
  }else{
    tmp$pi <- additive_logistic(tmp$eta)
  }

  # calc posterior probs and pi.
  if(nArchetypes>1){
    fits <- tmp$coefs
    logls_mus <- get_logls_sam(y = y, X = X, W = W, U = U, G = G,
                               S = S, spp_weights = spp_weights,
                               site_spp_weights = site_spp_weights,
                               offset = offset, y_is_na = y_is_na,
                               disty = disty, linky=linky, size = size,
                               powers = powers, control=control,
                               fits = fits, get_fitted = FALSE)
    tmp$tau <- get_taus(tmp$pi,logls_mus$logl_sp, G, S)
    tmp$pi <- colSums(tmp$tau)/S
  }

  # Information criteria
  tmp <- calc_info_crit_sam(tmp)

  # titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_sam(titbits, y, X, W, U, spp_weights,
                                 site_spp_weights, offset, y_is_na, size, powers,
                                 archetype_formula, species_formula, all_formula,
                                 data, control, family, linky)

  # remove large annoying object if titbits == FALSE
  if(!titbits)
    tmp <- titbit_cleanup(tmp)

  tmp2 <- sam_model_clean_up_names(tmp)
  tmp2$call <- call
  tmp2$terms <- tt

  class(tmp2) <- c("species_mix")
  return(tmp2)
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
#'@param linky The link function to use alongside the distribution.
#'@param size The size of each of binomial sample at each site. Length should be the number of sites.
#'@param powers The power parameters for the Tweedie distribution.
#'@param control this is a list of control parameters that alter the specifics of model fitting.
#'@param inits This will be a vector of starting values for species_mix (i.e you've fitted a model and want to refit it).
#'@export

"species_mix.fit" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
                              offset, y_is_na, disty, linky,
                              size, powers, control, inits=NULL){

  if(G==1){
    tmp <- fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights,
                       offset, y_is_na, G, S, disty, linky, size, powers,
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
                                                linky = linky,
                                                size = size,
                                                powers = powers,
                                                control = control)

  } else {
    if(!control$quiet)message('Be careful! You are using your own initial starting values to optimise the species_mix model.')
    inits <- setup_inits_sam(inits, S=S, G=G, X=X, W=W, U=U, disty, return_list = TRUE)
    # print(inits)
    starting_values <- inits
  }

  tmp <- sam_optimise(y = y, X = X, W = W, U = U, offset = offset,
                      spp_weights =  spp_weights,
                      site_spp_weights = site_spp_weights,
                      y_is_na = y_is_na, S = S, G = G, disty = disty, linky = linky,
                      size = size, powers = powers, start_vals = starting_values,
                      control = control)

  return(tmp)
}

#'@rdname species_mix.multifit
#'@name species_mix.multifit
#'@title species_mix.multifit
#' @details species_mix.multifit is used to fit species archetype models with
#' multiple fits. This can be useful for complex models where there are many
#' local maximum in the likelihood. The multiple restarts using in
#' `species_mix.multifit` can help search the log-likelihood using starting values
#' with some random variance (noise).
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
#'  you might use this an alternative to an offset, where there might be a set
#'  covariates which means artifact of how the data was collected or some other
#'  process like a seasonality.
#' @param data a matrix or data.frame which contains the 'species_data'
#' matrix, a const and the covariates in the structure of spp1, spp2, spp3,
#' const, temperature, rainfall. dims of matrix should be
#' nsites*(nspecies+const+covariates).
#' @param nArchetypes The number of archetypes (mixing components/groups) to
#' estimate from the data. This need to be explicitly declared. By default it is
#' set to 3.
#' @param family The family of statistical family to use within
#' the ecomix models. a  choice between "bernoulli", "binomial", "poisson",
#' "negative.binomial", "tweedie" and "gaussian" families are possible and
#' applicable to specific types of data. If you wish to fit a point process
#' version of a SAM please refer to the \link{ecomix.ppm} package.
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
#' @param power The power parameter for 'tweedie' family. Default is 1.6, and
#' this is assigned to all species
#' @param control a list of control parameters for optimization and calculation.
#' See details.
#' @param inits NULL a numeric vector that provides approximate starting values
#' for species_mix coefficients. These are family specific, but at a
#' minimum you will need pi (additive_logistic transformed), alpha
#' (intercepts) and beta (mixing coefs).
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
#'@description This version of species mix is useful for fitting models which
#'have complex likelihoods. The multiple starts will enable optimization of the
#'log-likelihood using multiple starts.
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
#' simulated_data <- species_mix.simulate(archetype_formula = sam_form,
#' species_formula = sp_form,
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
                                   control=list(ecm.prefit=FALSE),
                                   inits=NULL,
                                   titbits = FALSE,
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

  groupselection <- FALSE
  if(length(nArchetypes)>1){
    if(!control$quiet)
      message(paste0("Running species_mix.multifit to do groups selection from ",nArchetypes[1]," to ",nArchetypes[length(nArchetypes)], " archetypes\nUsing ",nstart," starts per archetype."))
      groupselection <- TRUE
   } else if (length(nArchetypes)==1){
     if(!control$quiet)
      message(paste0("Running species_mix.multifit using ",nstart," starts to fit ",nArchetypes," Archetypes to the data."))
   } else {
     stop("Please provide a single 'nArchetypes' value of a vector of integers to ranging from smallest to largest number of groups")
   }

  which_mix <- check_species_formula(species_formula)
  if(!is.null(species_formula))
    species_formula <- stats::as.formula(species_formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("data","offset","weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.exclude"
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # need this for the na.omit step
  rownames(mf)<-seq_len(nrow(mf))
  dat <- clean_data_sam(mf, archetype_formula, species_formula, all_formula, family)

  # get responses
  y <- stats::model.response(dat$mf.X)

  # logical matrix needed for removing NAs from response and weights.
  y_is_na <- is.na(y)

  # check names of reponses
  S <- check_reponse_sam(y)

  if (!S){
    if(!control$quiet)
      message("Two species have the same name -- exitting now")
    return(NULL)
  }
  if( !control$quiet)
    message( "There are ", paste(nArchetypes,collapse=" "), " archetypes")

  # get archetype model matrix
  Xdat <- get_X_sam(mf.X = dat$mf.X)
  X <- Xdat$X
  xterms <- Xdat$mt.x

  # what is the W matrix (species covariates)
  Wdat <- get_W_sam(mf.W = dat$mf.W)
  W <- Wdat$W
  wterms <- Wdat$mt.w

  # what is the U matrix (all covariates)
  Udat <- get_U_sam(mf.U = dat$mf.U)
  U <- Xdat$U
  uterms <- Xdat$mt.u

  tt <- list(xterms=xterms,wterms=wterms,uterms=uterms)

  # get offsets
  offset <- get_offset_sam(dat$mf.X)

  # get the weights
  species_names <- colnames(y)
  site_spp_weights <- get_site_spp_weights_sam(mf,weights,species_names)#,family)
  spp_weights <- check_spp_weights(bb_weights,S)

  # get size for binomial
  size <- check_size_binomial(size,nrow(dat$mf.X))

  # get family
  disty.listy <- get_family_sam(family,size)
  disty <- disty.listy$disty
  linky <- disty.listy$link
  family <- disty.listy$family

  # check powers
  powers <- get_power_sam(disty,power,S)

  # summarising data to console
  if(!control$quiet) print_input_sam(y, X, W, U, S, archetype_formula,
                                     species_formula, all_formula,
                                     family, linky,
                                     quiet=control$quiet)


  tmp_fun <- function(x,gg){
    if( !control$quiet & nstart>1)
      setTxtProgressBar(pb, x)
    # tmpQuiet <- control$quiet
    control$quiet <- TRUE;
    control$trace <- 0
    control$init.method <- "random2"
    # fit species mix.
    G <- nArchetypes <- gg
    tmp <- species_mix.fit(y=y, X=X, W=W, U=U, G=G, S=S,
                           spp_weights=spp_weights,
                           site_spp_weights=site_spp_weights,
                           offset=offset, disty=disty, linky=linky, y_is_na=y_is_na,
                           size=size, powers=powers, control=control, inits=inits)

    tmp$family <- disty_cases[disty]

    tmp$S <- S;
    tmp$G <- G;
    tmp$npx <- ncol(X);
    tmp$npw <- ifelse(ncol(W)>1,ncol(W),0);
    tmp$npu <- ifelse(!is.null(U),ncol(U),0);

    if(nArchetypes==1){
      tmp$pi <- tmp$pi
    }else{
      tmp$pi <- additive_logistic(tmp$eta)
    }

    #get logls from parameters
    #calc posterior probs and pi.
    if(nArchetypes>1){
      fits <- tmp$coefs
      logls_mus <- get_logls_sam(y = y, X = X, W = W, U = U, G = G,
                                 S = S, spp_weights = spp_weights,
                                 site_spp_weights = site_spp_weights,
                                 offset = offset, y_is_na = y_is_na,
                                 disty = disty, linky=linky, size = size,
                                 powers = powers, control=control,
                                 fits = fits, get_fitted = FALSE)
      tmp$tau <- get_taus(tmp$pi,logls_mus$logl_sp, G, S)
      tmp$pi <- colSums(tmp$tau)/S
    }

    #Information criteria
    tmp <- calc_info_crit_sam(tmp)

    #titbits object, if wanted/needed.
    tmp$titbits <- get_titbits_sam(titbits, y, X, W, U, spp_weights,
                                   site_spp_weights, offset, y_is_na, size, powers,
                                   archetype_formula, species_formula, all_formula,
                                   data, control, family)

    # remove large annoying object if titbits == FALSE
    if(!titbits)
      tmp <- titbit_cleanup(tmp)

    tmp2 <- sam_model_clean_up_names(tmp)
    tmp2$call <- call
    tmp2$terms <- tt

    class(tmp2) <- c("species_mix")
    return(tmp2)
  }


  if(groupselection){

    many_starts <- list()
    for(ii in nArchetypes){

      if( !control$quiet & nstart>1){
        message(paste0("\nFitting ",ii," Archetypes"))
        pb <- txtProgressBar(min = 1, max = nstart, style = 3)
      }

      many_starts[[ii]] <- plapply(seq_len(nstart), tmp_fun, gg = nArchetypes[ii],
                                .parallel = mc.cores,
                                .verbose = !control$quiet)

    }

  } else {
    if( !control$quiet & nstart>1)
        pb <- txtProgressBar(min = 1, max = nstart, style = 3)

   #Fit the model many times
   many_starts <- plapply(seq_len(nstart), tmp_fun, gg = nArchetypes,
                         .parallel = mc.cores,
                          .verbose = !control$quiet)

   }

   res <- list(multiple_fits = many_starts,
              nArchetypes=nArchetypes,
              groupselection=groupselection,
              nstart = nstart)


   class(res) <- c("species_mix.multifit")

   return(res)
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

  disty.listy <- get_family_sam(family,size)
  disty <- disty.listy$disty
  linky <- disty.listy$link
  family <- disty.listy$family

  link <- make.link(link = linky)


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

  # if(family %in% c('bernoulli','binomial')) link <- make.link('logit')
  # if(family %in% c('poisson','negative.binomial','tweedie')) link <- make.link('log')
  # if(family %in% c('gaussian')) link <- make.link('identity')
  # if(family %in% 'ippm') {
  #   grid <- simulate_ippm_grid(X,W)
  #   grid2D <- grid$grid2D
  #   X <- grid$X
  #   W <- grid$W
  # }

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
  attr(res, "pi") <- pi
  attr(res, "alpha") <- alpha
  attr(res, "beta") <- beta
  attr(res, "gamma") <- gamma
  attr(res, "delta") <- delta
  attr(res, "logTheta") <- logTheta
  attr(res, "powers") <- powers
  attr(res, "mu") <- fitted
  attr(res, "ippm_weights") <- wts
  attr(res, "size") <- size
  attr(res, "offset") <- offset
  return(res)
}


###### internal fitting functions #####
"fit.ecm.sam" <- function(y, X, W, U=NULL, spp_weights, site_spp_weights, offset,
                          y_is_na, G, S, disty, linky, size, powers, control, starting.sam = NULL){

  n <- nrow(y)
  starting.sam <-  NULL
  bestOfAllMods <- list( logl=-Inf)
  for(t in seq_len(control$ecm.refits)) {
    if(!control$quiet) message(paste0("ECM restart ",t," of ", control$ecm.refits,""))
    starting.sam <- get_initial_values_sam(y = y, X = X, W = W, U = U,
                                           site_spp_weights = site_spp_weights,
                                           offset = offset, y_is_na = y_is_na,
                                           G = G, S = S,
                                           disty = disty, linky = linky, size = size,
                                           control = control)


    fits <- list()
    tau <- starting.sam$tau
    fits$pi <- colMeans(tau)
    fits$alpha <- starting.sam$alpha
    fits$beta <- starting.sam$beta
    fits$gamma <- starting.sam$gamma
    fits$delta <- starting.sam$delta
    fits$theta <- starting.sam$theta

    diff.logl <- -100; cw.logl <- -Inf; new.logl <- 10
    mod <- list(logl = -Inf)

    counter <- 1
    while((diff.logl < 1-control$ecm.reltol | diff.logl > 1+control$ecm.reltol) & counter <= control$ecm.steps) {
      # if(diff.logl > 1) break;
      ## If any pi hits 0, restart with new random starts
    if(any(fits$pi < 0.0005)) {
    starting.sam <- get_initial_values_sam(y = y, X = X, W = W, U = U,
                                           site_spp_weights = site_spp_weights,
                                           offset = offset, y_is_na = y_is_na,
                                           G = G, S = S,
                                           disty = disty, linky = linky, size = size,
                                           control = set_control_sam(list(init.method="random2")))

      fits <- list()
      tau <- starting.sam$tau
      fits$pi <- colMeans(tau)
      fits$alpha <- starting.sam$alpha
      fits$beta <- starting.sam$beta
      fits$gamma <- starting.sam$gamma
      fits$delta <- starting.sam$delta
      fits$theta <- starting.sam$theta
      }

      ## CM-step (M step)
      fits$pi <- colMeans(tau)

      ## optimise the spp coefs - alpha & gamma if present
      fm_sppParam <- plapply(seq_len(S), apply_optimise_sppParam,
                             y, X, W, U, S, G, n, tau,
                             fits, site_spp_weights,
                             offset, y_is_na, disty, linky,
                             size, powers,
                             .parallel = control$cores,
                             .verbose = FALSE)
      new.sp.params <- do.call(rbind,lapply(fm_sppParam, `[[`, 1))
      fits$alpha <- update_coefs(fits$alpha,new.sp.params[,1])
      if(ncol(W)>1){
        fits$gamma <- update_coefs(fits$gamma,new.sp.params[,-1])
      } else {
        fits$gamma <- rep(-999999,S)
      }

      ## Update archetype-specific coefficients
      fm_beta <- plapply(seq_len(G), apply_optimise_betas,
                         y, X, W, U, site_spp_weights, offset,
                         y_is_na, disty, linky,
                         tau, fits, size, powers, control,
                         .parallel = control$cores,
                         .verbose = FALSE)
      new.beta <- do.call(rbind,fm_beta)
      fits$beta <- update_coefs(fits$beta,new.beta)

      ## Update the all parameter.
      if(!is.null(U)){
        new.delta <- stats::nlminb(start = fits$delta, objective = llogl.allParams,
                                   y = y, X = X, W = W, U=U, tau = tau, fits = fits,
                                   site_spp_weights = site_spp_weights,
                                   offset = offset, y_is_na = y_is_na,
                                   disty = disty, linky = linky,
                                   size = size, powers = powers)
        fits$delta <- update_coefs(fits$delta, new.delta)
      }

      ## update the dispersion parameters.
      if(disty %in% c(4,5,6)) {
        fm_theta <- plapply(seq_len(S), apply_optimise_thetaParams,
                                     y, X, W, U, tau, fits,
                                     site_spp_weights, offset, y_is_na,
                                     disty, linky, size, powers,
                                     .parallel = control$cores,
                                     .verbose = FALSE)
        new.sp.thetas <- unlist(lapply(fm_theta, `[[`, 1))
        fits$theta <- update_coefs(fits$theta,new.sp.thetas)
      } else {
        fits$theta <- rep(-999999,S)
      }

      ## E-step
      do.estep <- try(e.step(y, X, W, U, site_spp_weights,
                             offset, y_is_na, disty, linky,
                             fits, size, powers, control),
                      silent=TRUE)
      if(!is.finite(do.estep$logl)) {
        mod <- list(logl = -Inf); break;
      }
      if(is.finite(do.estep$logl)) {
        tau <- do.estep$tau
        new.logl <- do.estep$logl
        diff.logl <- abs(new.logl/cw.logl)

        if(!control$quiet) message(paste0("Iteration: ", counter, " | New loglik ", round(new.logl,3),  " | Ratio loglik ", round(diff.logl,6)));

        cw.logl <- new.logl
        counter <- counter + 1
      }
    }

    if(new.logl > mod$logl) {
      mod <- list(logl = new.logl, alpha = fits$alpha, beta = fits$beta,
                  eta = additive_logistic(fits$pi, inv = TRUE)[-G],
                  gamma = fits$gamma, delta = fits$delta,
                  theta = fits$theta, pi = fits$pi, tau = tau)
      mod$diff.logl <- diff.logl
      mod$counter <- counter
      names(mod$pi) <- paste("Archetype", 1:G, sep = "");
      mod$tau <- as.data.frame(round(mod$tau,5)); colnames(mod$tau) <- names(mod$pi); rownames(mod$tau) <- paste("Sp",1:S,sep="")
      mod$entropy <- -sum(as.vector(mod$tau[mod$tau != 0])*log(unlist(mod$tau[mod$tau != 0])));
      mod$beta <- round(mod$beta,6); rownames(mod$beta) <- names(mod$pi);
      names(mod$alpha) <- rownames(mod$tau)
      if(ncol(W)>1) {
        rownames(gamma) <- rownames(mod$tau)
        colnames(gamma) <- colnames(W[,-1,drop=FALSE])
      } else {
        names(mod$gamma) <- rownames(mod$tau)
      }
      if(!is.null(U)) names(mod$beta) <- colnames(U)
      names(mod$theta) <- rownames(mod$tau)
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
                     y_is_na, disty, linky,
                     fits, size, sp.powers, control,
                     get.fitted = FALSE) {
  S <- ncol(y); n <- nrow(y); G <- length(fits$pi)
  out.tau <- matrix(0,S,G)
  sp.logl <- rep(0,S) ## Species specific incomplete logL

  link <- make.link(link = linky)
  # if(disty%in%c(1,7)) link <- make.link('logit')
  # if(disty%in%c(2,3,4,5)) link <- make.link('log')
  # if(disty%in%6) link <- make.link('identity')

  if(get.fitted) fitted.values <- array(0,dim=c(G,n,S))

  if(!is.null(U)){
    all.etas <- U%*%c(fits$delta)
  } else {
    all.etas <- rep(0,n)
  }

  for(gg in seq_len(G)) {
    mix.etas <- X%*%fits$beta[gg,]

    for(ss in seq_len(S)) {
      sp_idx <- !y_is_na[,ss]
      if(ncol(W)>1){
        spp.etas <- W %*% c(fits$alpha[ss],fits$gamma[ss,])
      }else{
        spp.etas <- W %*% c(fits$alpha[ss])
      }
      new.etas <- all.etas + mix.etas + spp.etas + offset
      if(disty %in% 1) {
        check.p <- link$linkinv(new.etas)
        check.p[check.p < 1e-4] <- 1e-4; check.p[check.p > (1-1e-4)] <- (1-1e-4);
        if(get.fitted) fitted.values[gg,,ss] <- check.p
        out.tau[ss,gg] <- (sum(dbinom(y[,ss], 1, prob = check.p, log = TRUE)))
      }
      if(disty %in% 2) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.tau[ss,gg] <- (sum(dpois(y[,ss], lambda = link$linkinv(new.etas), log = TRUE)))
      }
      # if(disty %in% 3) {
      #   if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
      #   out.tau[ss,gg] <- (y[sp_idx,ss] %*% new.etas[sp_idx,,drop=FALSE] - site_spp_weights[sp_idx,ss] %*% link$linkinv(new.etas[sp_idx,,drop=FALSE]))
      # }
      if(disty %in% 4) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.tau[ss,gg] <- (sum(dnbinom(y[,ss], mu = link$linkinv(new.etas), size = 1/fits$theta[ss], log = TRUE)))
      }
      if(disty %in% 5) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.tau[ss,gg] <- (sum(fishMod::dTweedie(y[,ss], mu = link$linkinv(new.etas), phi = fits$theta[ss], p = sp.powers[ss], LOG = TRUE)))
      }
      if(disty %in% 6) {
        if(get.fitted) fitted.values[gg,,ss] <- link$linkinv(new.etas)
        out.tau[ss,gg] <- (sum(dnorm(y[,ss], mean = link$linkinv(new.etas), sd = sqrt(fits$theta[ss]), log = TRUE)))
      }
      if(disty %in% 7) {
        check.p <- link$linkinv(new.etas)
        check.p[check.p < 1e-4] <- 1e-4; check.p[check.p > (1-1e-4)] <- (1-1e-4);
        if(get.fitted) fitted.values[gg,,ss] <- check.p
        out.tau[ss,gg] <- (sum(dbinom(y[,ss], size, prob = check.p, log = TRUE)))
      }

    }
  }

  for(ss in seq_len(S)) {
    eps <- max(out.tau[ss,])
    sp.logl[ss] <- log(sum(fits$pi*exp(out.tau[ss,]-eps))) + eps
  }
  for(ss in seq_len(S)) {
    for(gg in seq_len(G)) {
      out.tau[ss,gg] <- exp((log(fits$pi[gg]) + out.tau[ss,gg]) - sp.logl[ss])
    }
  }

  full.logl <- sum(sp.logl)

  ## Trying to add in the penalties.
  if(control$doPenalties>0){
    alphaPen <- calc_alpha_penalties(fits,control)
    betaPen <- calc_beta_penalties(fits,control)
    if(ncol(W)>1){
      gammaPen <- calc_gamma_penalties(fits,control)
    } else {
      gammaPen <- 0
    }
    if(!is.null(U)){
      deltaPen <- calc_delta_penalties(fits,control)
    } else {
      deltaPen <- 0
    }
    if(disty%in%c(4,5,6)){
      thetaPen <- calc_theta_penalties(fits,control)
    } else {
      thetaPen <- 0
    }

    penalties <- alphaPen + betaPen + gammaPen + deltaPen + thetaPen

    full.logl <- full.logl + penalties

  }

  out.list <- list(tau = out.tau, sp.logl = sp.logl, logl = full.logl)
  if(get.fitted) out.list$fitted = fitted.values
  return(out.list)
}

"apply_optimise_sppParam" <- function(ss, y, X, W, U, S, G, n, tau, fits,
                                      site_spp_weights, offset, y_is_na,
                                      disty, linky, size, powers){

  if(ncol(W)>1) pars <- c(fits$alpha[ss],fits$gamma[ss,])
  else pars <- fits$alpha[ss]

  update.sppParams <- stats::nlminb(start = pars, objective = llogl.sppParams, ss = ss,
                               y = y, X = X, W = W, U=U, S=S, G=G, n=n, tau = tau, fits = fits,
                               site_spp_weights = site_spp_weights,
                               offset = offset, y_is_na = y_is_na,
                               disty = disty, linky = linky, size = size, powers = powers)
  return(update.sppParams)
}

"apply_optimise_thetaParams" <- function(ss, y, X, W, U, tau, fits,
                                   site_spp_weights, offy, y_is_na,
                                   disty, linky, size, powers){
  update.theta <- suppressWarnings(stats::optimize(f = llogl.thetaParams, interval = c(0.001,100),
                                                   maximum = TRUE, ss = ss,
                           y = y, X = X, W = W, U=U, tau = tau, fits = fits,
                           site_spp_weights = site_spp_weights,
                           offset = offy, disty = disty, linky=linky,
                           size = size, powers = powers)$maximum)
  return(update.theta)
}

"apply_optimise_betas" <- function(gg, y, X, W, U, site_spp_weights, offset,
                                   y_is_na, disty, linky, tau,
                                   fits, size, powers, control){

  n <- nrow(y)
  if(disty %in% c(2,3,4,5)){
  glmnet.family <- "poisson"; glm.family <- poisson();
  }
  if(disty %in% c(1,7)){
    if(linky=="logit")
      glmnet.family <- glm.family <- "binomial"
    if(linky=="cloglog")
      glmnet.family <- glm.family <- binomial(link="cloglog");
  }
  if(disty %in% c(6)){
  glmnet.family <- "gaussian"; glm.family <- gaussian();
  }

  if(disty %in% 4) get.mus <- e.step(y, X, W, U, site_spp_weights,
                                     offset, y_is_na, disty, linky,
                                     fits, size, powers, control,get.fitted = TRUE)$fitted

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
  tau.weights <- rep(tau[,gg,drop=FALSE],c(n_ys))
        if(disty %in% 4) tau.weights <- rep(tau[,gg,drop=FALSE],c(n_ys))/(1+rep(fits$theta,each=n)*as.vector(get.mus[gg,,]))
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

        if(disty%in%c(7)){
          Y_s <- as.matrix(cbind(Y_s,size_s-Y_s))
        }

        if(ncol(X_s)==1) X_s <- cbind(1,X_s)

        if (disty%in%c(7)){
          fit1 <- try(glm2::glm.fit2(y=Y_s, x=as.data.frame(cbind(1,X_s)),weights=c(obs.weights),
                               family=glm.family, offset=offy), silent=FALSE)
          if (any(class(fit1)[1] %in% 'try-error')){
            new.betas <- rep(NA, ncol(X))
            names(new.betas) <- colnames(X)
          } else {
             new.betas <- coef(fit1)[-1]
          }
        } else {
          fit1 <- glmnet::glmnet(x = as.matrix(X_s), y = Y_s, family = glmnet.family, weights = obs.weights+1e-9, offset = offy, nlambda = 100, intercept = FALSE)
          new.betas <- coef(fit1)[,ncol(coef(fit1))][-1]
        }
        # print(new.betas)
        if(ncol(X)==1)new.betas <- new.betas[-1]

        return(new.betas)

}

"calc_alpha_penalties" <- function(fits, control){

  alphaPen <- -1 * fits$alpha*fits$alpha/(2*control$penalty.alpha*control$penalty.alpha)
  alphaPenSum <- sum(alphaPen)
  return(alphaPenSum)

}

"calc_beta_penalties" <- function(fits, control){

  betaPen <- -1 * fits$beta*fits$beta/(2*control$penalty.beta*control$penalty.beta)
  betaPenSum <- sum(betaPen)
  return(betaPenSum)

}

"calc_gamma_penalties" <- function(fits, control){

  gammaPen <- -1 * fits$gamma*fits$gamma/(2*control$penalty.gamma*control$penalty.gamma)
  gammaPenSum <- sum(gammaPen)
  return(gammaPenSum)

}

"calc_delta_penalties" <- function(fits, control){

  deltaPen <- -1 * fits$delta*fits$delta/(2*control$penalty.delta*control$penalty.delta)
  deltaPenSum <- sum(deltaPen)
  return(deltaPenSum)

}

"calc_theta_penalties" <- function(fits, control){

  thetaPen <- -1 * (log(1/fits$theta)*control$penalty.theta[1])*(log(1/fits$theta)*control$penalty.theta[1])/(2*control$penalty.theta[2]*control$penalty.theta[2])
  thetaPenSum <- sum(thetaPen)
  return(thetaPenSum)

}


"llogl.sppParams" <- function(x, ss, y, X, W, U, S, G, n, tau,
                              fits, site_spp_weights, offset,
                              y_is_na, disty, linky, size, powers) {

  # S <- ncol(y); n <- nrow(y); G <- ncol(tau)
  out <- 0
  sp_idx <- !y_is_na[,ss]

  link <- make.link(link=linky)
  # if(disty%in%c(1,7)) link <- make.link('logit')
  # if(disty%in%c(2,3,4,5)) link <- make.link('log')
  # if(disty%in%6) link <- make.link('identity')

  if(!is.null(U)) all.etas <- U%*%c(fits$delta)
  else all.etas <- rep(0,n)

  mix.etas <- tcrossprod(X,fits$beta)

  for(gg in seq_len(G)) {
    eta <- W%*%x + mix.etas[,gg] + all.etas + offset
    if(disty %in%  1)
      out <- out + sum(tau[ss,gg]*dbinom(y[,ss], size = 1, prob =  link$linkinv(eta), log = TRUE))
    if(disty %in%  2)
      out <- out + sum(tau[ss,gg]*dpois(y[,ss], lambda =  link$linkinv(eta), log = TRUE))
    # if(disty %in%  3)
    #   out <- out + sum(tau[ss,gg]*(dpois(y[sp_idx,ss]/site_spp_weights[sp_idx,ss], lambda =  link$linkinv(eta), log = TRUE)*site_spp_weights[sp_idx,ss]))#(y[sp_idx,ss] %*% eta[sp_idx,,drop=FALSE] - site_spp_weights[sp_idx,ss] %*% link$linkinv(eta[sp_idx,,drop=FALSE])))
    if(disty %in%  4)
      out <- out + sum(tau[ss,gg]*dnbinom(y[,ss], mu =  link$linkinv(eta), size = 1/fits$theta[ss], log = TRUE))
    if(disty %in%  5)
      out <- out + sum(tau[ss,gg]*fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
    if(disty %in%  6)
      out <- out + sum(tau[ss,gg]*dnorm(y[,ss], mean =  link$linkinv(eta), sd = sqrt(fits$theta[ss]), log = TRUE))
    if(disty %in%  7)
      out <- out + sum(tau[ss,gg]*dbinom(y[,ss], size = size, prob = link$linkinv(eta), log = TRUE))
  }
  return(-out) #negative for nlmnib
}


"llogl.thetaParams" <- function(x, ss, y, X, W, U, tau, fits,
                                site_spp_weights, offset,
                                disty, linky, size, powers) {

  S <- ncol(y); n <- nrow(y); G <- ncol(tau)
  out <- 0
  link <- make.link(link = linky)
  # if(disty%in%c(4,5)) link <- make.link('log')
  # if(disty%in%6) link <- make.link('identity')


  if(!is.null(U)) all.etas <- U%*%c(fits$delta)
  else all.etas <- rep(0,n)

  if(ncol(W)>1){
    spp.etas <- W %*% c(fits$alpha[ss],fits$gamma[ss,])
  }else{
    spp.etas <- W %*% c(fits$alpha[ss])
  }

  mix.etas <- tcrossprod(X,fits$beta)

  for(gg in seq_len(G)){
    eta <- spp.etas + mix.etas[,gg] + all.etas + offset
   if(disty==4)
     cw.out <- sum(tau[ss,gg]*dnbinom(y[,ss], mu = link$linkinv(eta), size = 1/x, log = TRUE));
   if(disty==5)
     cw.out <- sum(tau[ss,gg]*fishMod::dTweedie(y[,ss], mu = link$linkinv(eta), phi = x, prob = powers[ss], LOG = TRUE));
   if(disty==6)
     cw.out <- sum(tau[ss,gg]*dnorm(y[,ss], mean = eta, sd = sqrt(x), log = TRUE))
   out <- out + cw.out

  }
  return(out)
}



"llogl.allParams" <- function(x, y, X, W, U, tau, fits,
                              site_spp_weights, offset, y_is_na,
                              disty, linky, size, powers) {

  S <- ncol(y); n <- nrow(y); G <- ncol(tau)
  out <- 0

  link <- make.link(link = linky)
  # if(disty%in%c(1,7)) link <- make.link('logit')
  # if(disty%in%c(2,3,4,5)) link <- make.link('log')
  # if(disty%in%6) link <- make.link('identity')

  all.etas <- U%*%x
  mix.etas <- tcrossprod(X,fits$beta)

  for (ss in seq_len(S)){

    sp_idx <- !y_is_na[,ss]
    if(ncol(W)>1){
      spp.etas <- W %*% c(fits$alpha[ss],fits$gamma[ss,])
    }else{
      spp.etas <- W %*% c(fits$alpha[ss])
    }

    for(gg in seq_len(G)) {
      eta <- spp.etas + mix.etas[,gg] + all.etas + offset
      if(disty %in%  1)
        out <- out + sum(tau[ss,gg]*dbinom(y[,ss], size = 1, prob =  link$linkinv(eta), log = TRUE))
      if(disty %in%  2)
        out <- out + sum(tau[ss,gg]*dpois(y[,ss], lambda =  link$linkinv(eta), log = TRUE))
      # if(disty %in%  3)
      #   out <- out + sum(tau[ss,gg]*(y[sp_idx,ss] %*% eta[sp_idx,,drop=FALSE] - site_spp_weights[sp_idx,ss] %*%  link$linkinv(eta[sp_idx,,drop=FALSE])))
      if(disty %in%  4)
        out <- out + sum(tau[ss,gg]*dnbinom(y[,ss], mu =  link$linkinv(eta), size = 1/fits$theta[ss], log = TRUE))
      if(disty %in%  5)
        out <- out + sum(tau[ss,gg]*fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], prob = powers[ss], LOG = TRUE))
      if(disty %in%  6)
        out <- out + sum(tau[ss,gg]*dnorm(y[,ss], mean =  link$linkinv(eta), sd = sqrt(fits$theta[ss]), log = TRUE))
      if(disty %in%  7)
        out <- out + sum(tau[ss,gg]*dbinom(y[,ss], size = size, prob = link$linkinv(eta), log = TRUE))
    }

  }
  return(out)
}


"starting_values_wrapper" <- function(y, X, W, U, spp_weights, site_spp_weights,
                                      offset, y_is_na, G, S, disty, linky, size, powers, control){
  if(control$ecm.prefit){
    if(!control$quiet)message('Using ECM algorithm to find starting values; using ',
                              control$ecm.refits,' refits')
    emfits <- fit.ecm.sam(y, X, W, U, spp_weights, site_spp_weights,
                          offset, y_is_na, G, S, disty, linky, size,
                          powers, control)
    start_vals <- list(alpha = (emfits$alpha),
                       beta = (emfits$beta),
                       gamma = (emfits$gamma),
                       delta = (emfits$delta),
                       theta = (emfits$theta),
                       pi = (emfits$pi))
  } else {
    if(!control$quiet)message('If you are choosing to not use the ECM algorithm to estimate starting values\nwe recommend that you at least run a multifits to optimise the loglikelihood, see "species_mix.multifit".')
    if(!control$quiet)message('You are not using the ECM algorith to find starting values;\n starting values are generated using ',control$init.method,'.')
    starting_values <- get_initial_values_sam(y = y, X = X, W= W, U = U,
                                              site_spp_weights = site_spp_weights,
                                              offset = offset, y_is_na = y_is_na,
                                              G = G, S = S,
                                              disty=disty, linky=linky,
                                              size=size, powers = powers,
                                              control = control)
    start_vals <- list(alpha=(starting_values$alpha),
                       beta=(starting_values$beta),
                       gamma=(starting_values$gamma),
                       delta=(starting_values$delta),
                       theta=(starting_values$theta),
                       pi=(starting_values$pi))
  }

  ## all the things we need to c++ optimisation.
  start_vals$eta <- additive_logistic(start_vals$pi, inv = TRUE)[-G]
  start_vals$nS <- S
  start_vals$nG <- G
  start_vals$nObs <- nrow(y)
  #put dispersion parameters on the log-scale for c++
  if(disty%in%c(4,5,6)) start_vals$theta <- transform.theta(start_vals$theta, disty=disty, logify = TRUE)
  return(start_vals)
}

## function for starting values using penalities
"get_initial_values_sam" <- function(y, X, W, U=NULL, site_spp_weights,
                                   offset, y_is_na, G, S, disty, linky,
                                   size, powers, control) {

  n <- nrow(y)
  prev.min <- floor(n*control$minimum.sites);
  sel.omit.spp <- which(colSums(y>0) <= prev.min)
  if(length(sel.omit.spp)==0) sel.omit.spp <- -1*(1:S)

  starting.sam <- list(alpha = rep(0,S), theta = rep(1,S));
  if(!control$quiet) message("Initialising starting values")


  if(is.null(U)){
    datX <- cbind(W,X)
  } else {
    datX <- cbind(W,X,U)
  }

  datX <- as.data.frame(datX)

  tmpform <- as.formula(paste0("mvabund::mvabund(y) ~ - 1 + ."))

  if(disty%in%1){ # bernoulli
    # if(linky==)
    fam <- binomial(link = linky)
    # if(linky==1) fam <- binomial(link = 'cloglog')
    fit1 <- mvabund::manyglm(tmpform, data = datX, offset = offset, family = fam)
    coefs <- t(fit1$coefficients)
    starting.sam$theta <- rep(-999999,S)
  }
  if(disty%in%2){ # poisson
    fit1 <- mvabund::manyglm(tmpform, data = datX, offset = offset, family = "poisson")
    coefs <- t(fit1$coefficients)
    starting.sam$theta <- rep(-999999,S)
  }
  # if(disty%in%3){ # ippm
  #   fit1 <- many.fit(y, X, W, U, site_spp_weights,
  #                     offset, y_is_na, G, S, disty, size, powers, control)
  #   coefs <- fit1$coefficients
  # }
  if(disty%in%4){ # negative binomial
    fit1 <- mvabund::manyglm(tmpform, data = datX,offset = offset, family = "negative.binomial")
    coefs <- t(fit1$coefficients)
    starting.sam$theta <- fit1$phi
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
    starting.sam$theta <- rep(-999999,S)
  }

  # a bit of book keeping
  covarNames <- colnames(datX)
  alphaIdx <- 1#match(colnames(W[,1,drop=FALSE]),covarNames)
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
    starting.sam$gamma <- -999999
  }

  # delta for entire dataset
  if(!is.null(U)){
    deltaIdx <- match(colnames(U),covarNames)
    starting.sam$delta <- colMeans(coefs[-sel.omit.spp,deltaIdx,drop=FALSE])
  } else {
    starting.sam$delta <- -999999
  }

  if(G==1) control$init.method <- 'kmeans'

  if(control$init.method=='kmeans'){
    if(!control$quiet) message("Initial groups parameter estimates by K-means clustering")
    fmmvnorm <- stats::kmeans(spp.beta, centers=G, iter.max = 200, nstart = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(spp.beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G)
    colnames(grp_coefs) <- covarNames[betaIdx]
    rownames(grp_coefs) <-  paste("Archetype", 1:G, sep = "");
  }

  if(control$init.method=='kmed'){
    if(!control$quiet) message("Initial groups parameter estimates by K-medoids")
    mrwdist <- kmed::distNumeric(spp.beta, spp.beta, method = "mrw")
    fmmvnorm <- kmed::fastkmed(mrwdist, ncluster = G, iterate = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- spp.beta[fmmvnorm$medoid,,drop=FALSE]
    grp_coefs <- matrix(grp_coefs,nrow=G)
    colnames(grp_coefs) <- covarNames[betaIdx]
    rownames(grp_coefs) <-  paste("Archetype", 1:G, sep = "");
  }

  if(control$init.method=='random2'){
    if(!control$quiet) message("Initial groups parameter estimates by K-means clustering with random noise")
    fmmvnorm <- stats::kmeans(spp.beta, centers=G, iter.max = 200, nstart = 100)
    tmp_grp <- fmmvnorm$cluster
    grp_coefs <- apply(spp.beta, 2, function(x) tapply(x, tmp_grp, mean))
    grp_coefs <- matrix(grp_coefs,nrow=G) # not check this...
    colnames(grp_coefs) <- covarNames[betaIdx]
    rownames(grp_coefs) <-  paste("Archetype", 1:G, sep = "");
    random_coefs <- sam_random_inits(alpha = starting.sam$alpha,beta = grp_coefs,
                                              gamma = starting.sam$gamma, delta = starting.sam$delta,
                                              theta = starting.sam$theta, S, G, X, W, U,
                                              disty, mult=0.3, control$init.sd)
    starting.sam$alpha <- random_coefs[[1]]
    grp_coefs <- random_coefs[[2]]
    starting.sam$gamma <- random_coefs[[3]]
    starting.sam$delta <- random_coefs[[4]]
    if(disty%in%c(4,5,6)) starting.sam$theta <- random_coefs[[5]]
  }

  #get tau as starting values
  if(G==1){
    tau <- matrix(1,nrow=ncol(y), ncol = G)
  } else {
    tau <- matrix(0,nrow=ncol(y), ncol = G)
    if(length(sel.omit.spp)>0){
      for(j in 1:length((1:S)[-sel.omit.spp]))
        tau[(1:S)[-sel.omit.spp][j],fmmvnorm$cluster[j]] <- 1
      tau[sel.omit.spp,] <- matrix(runif(length(sel.omit.spp)*G),length(sel.omit.spp), G)
    } else {
      for(j in seq_len(S))tau[j,fmmvnorm$cluster[j]] <- 1
    }
  }

  tau <- tau/rowSums(tau)
  # starting.sam$alpha
  starting.sam$beta <- grp_coefs
  if(!disty%in%3){
    starting.sam$tau <- shrink_taus(tau, G=G)
  } else {
    starting.sam$tau <- tau
  }
  starting.sam$pi <- colMeans(tau)

  return(starting.sam)
}

"many.fit" <- function(y, X, W, U, site_spp_weights, offset, y_is_na, G, S,
                       disty, linky, size, powers, control){

  options(warn = -1)
  fm_sp_mods <-  plapply(seq_len(S), apply_species_fits, y, X, W, U,
                         site_spp_weights, offset, y_is_na, disty, linky,
                         size, powers,
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
                                 offset, y_is_na, disty, linky, size, power){

  # which family to use?
  if(disty %in% c(1,7)) #binomials
    fam <- binomial(link=linky) #glmnet
  if(disty %in% c(2,3,4))
    fam <- "poisson"
  if(disty %in% 6)
    fam <- "gaussian"

  ids_i <- !y_is_na[,ss]

  # if (disty==3){
  #   outcomes <- as.numeric(y[ids_i,ss]/site_spp_weights[ids_i,ss])
  # } else {
    outcomes <- as.matrix(y[ids_i,ss])
  # }

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
    # if(disty==3) lambda.seq <- 0
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
      # if(tmp>5) tmp <- 5
      theta <- tmp
    }
    if( disty == 6){
      preds <- as.numeric( predict(ft_sp, s=locat.s, type="link",
                                   newx=as.matrix(df), newoffset=offset[ids_i]))
      #should be something like the resid standard
      theta <- sum((outcomes - preds)^2)/length(outcomes)
    }
  }

  if (disty==5) { #Tweedie needs an unconstrained fit.  May cause problems in some cases, especially if there is quasi-separation...
      df3 <- data.frame(y=outcomes, offy=offset, Intercept= 1, df)
      tmp.fm1 <- fishMod::tglm(y~-1+.-offy+offset( offy),
                               wts=c(site_spp_weights[,ss]),
                               data=df3, p=power[ss], vcov=FALSE,
                               residuals=FALSE, trace=0)
      my_coefs <- t(as.matrix(tmp.fm1$coef,ncol=1))
      theta <- tmp.fm1$coef["phi"]
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
                           S, G, disty, linky, size, powers, start_vals, control){

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
    gamma.score <- as.numeric(rep(NA, length(gamma)))
    Wcpp <- W[,-1,drop=FALSE]
  } else {
    npw <- as.integer(1)
    control$optiPart <- as.integer(0)
    gamma.score <- as.numeric(rep(NA, length(gamma)))
    Wcpp <- matrix(1,nrow = n, ncol=1)
  }

if(!is.null(U)) {
    Ucpp <- U
    npu <- as.integer(ncol(U))
    control$optiAll <- as.integer(1)
    delta.score <- as.numeric(rep(NA, length(delta)))
  } else {
    delta.score <- as.numeric(rep(NA, length(delta)))
    control$optiAll <- as.integer(0)
    Ucpp <- matrix(1,nrow = n,ncol=1)
    npu <- as.integer(1) # a dummy variable to stop c++ issues.
  }

  if(disty%in%c(4,5,6)){
    control$optiDisp <- as.integer(0)
    theta.score <- as.numeric(rep(NA, length(theta)))
  }else{
    control$optiDisp <- as.integer(0)
    theta.score <- as.numeric(rep(NA, length(theta)))
  }
  scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,
                                       gamma.score,delta.score,theta.score))))
  control$conv <- as.integer(0)

  if(linky=="cloglog") linkyin <- as.integer(1)
  else linkyin <- as.integer(0)

  #model quantities
  pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
  mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
  loglikeS <- as.numeric(rep(NA, S))
  loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))

  if(control$print.cpp.vals)print_starting_values(as.integer(S),
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
               as.integer(disty),as.integer(linkyin),
               as.integer(control$optiDisp),
               as.integer(control$optiPart),as.integer(control$optiAll),
               as.integer(control$doPenalties), ## should you do penalties?
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
               as.integer(control$optimise.cpp), as.integer(control$loglOnly.cpp),
               as.integer(control$derivOnly.cpp),
               # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
               PACKAGE = "ecomix")

  ret <- tmp
  ret$logl <- ret$logl * -1
  ret$mus <- array(mus, dim=c(n, S, G))
  ret$names <- list(spp=colnames(y), SAMs=paste("Archetype", seq_len(G), sep=""),
                    Xvars=colnames(X), Wvars=colnames(Wcpp), Uvars = colnames(Ucpp))

  if(disty%in%c(4,5,6)) ret$theta <- transform.theta(ret$theta,disty,logify=FALSE)

  names(ret$alpha) <- ret$names$spp
  names(ret$beta) <- paste(rep(ret$names$Xvars,each=G),ret$names$SAMs,sep='.')
  names(ret$eta) <- paste0("eta",seq_len(G-1))
  if(ncol(W)>1) names(ret$gamma) <- paste(rep(ret$names$Wvars,each=S),ret$names$spp,sep='.')
  if(!is.null(U)) names(ret$delta) <- ret$names$Uvars
  names(ret$theta) <- paste0("theta.",ret$names$spp)

  ret$coefs <- list(alpha = ret$alpha,
                    beta = matrix(ret$beta,G,npx),
                    eta = ret$eta,
                    gamma = matrix(ret$gamma,S,npw),
                    delta = ret$delta,
                    theta = ret$theta)

  ret$scores <- list(alpha.scores = alpha.score,
                     beta.scores = beta.score,
                     eta.scores=eta.score,
                     gamma.scores = gamma.score,
                     delta.scores=delta.score,
                     theta.scores=theta.score)

  ret$S <- S; ret$G <- G; ret$npx <- npx; ret$npw <- ifelse(ncol(W)>1,ncol(W),0);
  ret$npu <- ifelse(!is.null(U),ncol(U),0); ret$n <- n; ret$disty <- disty;
  ret$start.vals <- inits
  ret$loglikeSG <- matrix(loglikeSG,  nrow = S, ncol = G, byrow = FALSE) # loglikes got mixed up. Need to address in c++ code.  #for residuals
  ret$loglikeS <- loglikeS  #for residuals
  gc()
  return(ret)
}

# "sam_optimise_tweedie" <- function(y, X, W, U, offset, spp_weights, site_spp_weights, y_is_na,
#                                    S, G, disty, size, powers, start_vals, control){
#
#     Tw.phi.func <- function( phi1, spp3){
#       disp3 <- theta
#       disp3[spp3] <- phi1
#       tmp1 <- .Call("species_mix_cpp",
#                    as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
#                    as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
#                    as.numeric(offset), as.numeric(spp_weights),
#                    as.numeric(as.matrix(site_spp_weights)),
#                    as.integer(as.matrix(!y_is_na)),
#                    as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
#                    as.integer(npw), as.integer(npu), as.integer(n),
#                    as.integer(disty),as.integer(TRUE),
#                    as.integer(control$optiPart),as.integer(control$optiAll),
#                    # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
#                    as.double(alpha), as.double(beta), as.double(eta),
#                    as.double(gamma), as.double(delta), as.double(disp3),
#                    as.double(powers),
#                    # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
#                    as.numeric(control$penalty.alpha),as.numeric(control$penalty.beta),
#                    as.numeric(control$penalty.pi),as.numeric(control$penalty.gamma),
#                    as.numeric(control$penalty.delta),
#                    as.numeric(control$penalty.theta[1]),
#                    as.numeric(control$penalty.theta[2]),
#                    # SEXP &RalphaPen, SEXP &RbetaPen, SEXP &RpiPen,  SEXP &RgammaPen,
#                    # SEXP &RdeltaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen,
#                    alpha.score, beta.score, eta.score, gamma.score, delta.score,
#                    theta.score, as.integer(control$getscores.cpp), as.numeric(scores),
#                    # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
#                    pis_out, mus, loglikeS, loglikeSG,
#                    # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
#                    as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
#                    as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
#                    as.integer(control$printparams.cpp),
#                    # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
#                    as.integer(FALSE), as.integer(TRUE), as.integer(FALSE),
#                    # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
#                    PACKAGE = "ecomix")
#       return( -as.numeric( tmp1))
#     }
#
#     Tw.phi.func.grad <- function( phi1, spp3){
#       disp3 <- theta
#       disp3[spp3] <- phi1
#       tmp.disp.score <- rep( -99999, S)
#
#       tmp1 <- .Call("species_mix_cpp",
#                     as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
#                     as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
#                     as.numeric(offset), as.numeric(spp_weights),
#                     as.numeric(as.matrix(site_spp_weights)),
#                     as.integer(as.matrix(!y_is_na)),
#                     as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
#                     as.integer(npw), as.integer(npu), as.integer(n),
#                     as.integer(disty),as.integer(TRUE),
#                     as.integer(control$optiPart),as.integer(control$optiAll),
#                     # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
#                     as.double(alpha), as.double(beta), as.double(eta),
#                     as.double(gamma), as.double(delta), as.double(disp3),
#                     as.double(powers),
#                     # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
#                     as.numeric(control$penalty.alpha),as.numeric(control$penalty.beta),
#                     as.numeric(control$penalty.pi),as.numeric(control$penalty.gamma),
#                     as.numeric(control$penalty.delta),
#                     as.numeric(control$penalty.theta[1]),
#                     as.numeric(control$penalty.theta[2]),
#                     # SEXP &RalphaPen, SEXP &RbetaPen, SEXP &RpiPen,  SEXP &RgammaPen,
#                     # SEXP &RdeltaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen,
#                     alpha.score, beta.score, eta.score, gamma.score, delta.score,
#                     theta.score, as.integer(control$getscores.cpp), as.numeric(scores),
#                     # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
#                     pis_out, mus, loglikeS, loglikeSG,
#                     # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
#                     as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
#                     as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
#                     as.integer(control$printparams.cpp),
#                     # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
#                     as.integer(FALSE), as.integer(FALSE), as.integer(TRUE),
#                     # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
#                     PACKAGE = "ecomix")
#
#       return( -as.numeric( tmp.disp.score[spp3]))
#     }
#
#     inits <- unname(c(start_vals$alpha, start_vals$beta, start_vals$eta, start_vals$gamma,
#                       start_vals$delta, start_vals$theta))
#     npx <- as.integer(ncol(X))
#     n <- as.integer(nrow(X))
#
#     # parameters to optimise
#     alpha <- as.numeric(start_vals$alpha)
#     beta <- as.numeric(start_vals$beta)
#     eta <- as.numeric(start_vals$eta)
#     gamma <- as.numeric(start_vals$gamma)
#     delta <- as.numeric(start_vals$delta)
#     theta <- as.numeric(start_vals$theta)
#
#     #scores
#     getscores <- 1
#     alpha.score <- as.numeric(rep(NA, length(alpha)))
#     beta.score <- as.numeric(rep(NA, length(beta)))
#     eta.score <- as.numeric(rep(NA, length(eta)))
#
#     # sort of the species formula structures for cpp
#     if(ncol(W)>1){
#       npw <- as.integer(ncol(W[,-1,drop=FALSE]))
#       control$optiPart <- as.integer(1)
#       gamma.score <- as.numeric(matrix(NA, nrow=S, ncol=ncol(W)))
#       Wcpp <- W[,-1,drop=FALSE]
#     } else {
#       npw <- as.integer(1)
#       control$optiPart <- as.integer(0)
#       gamma.score <- -99999
#       Wcpp <- matrix(1,nrow = n, ncol=1)
#     }
#
#     if(!is.null(U)) {
#       Ucpp <- U
#       npu <- as.integer(ncol(U))
#       control$optiAll <- as.integer(1)
#       delta.score <- as.numeric(matrix(NA, ncol=ncol(U)))
#     } else {
#       delta.score <- -99999
#       control$optiAll <- as.integer(0)
#       Ucpp <- matrix(1,nrow = n,ncol=1)
#       npu <- as.integer(1) # a dummy variable to stop c++ issues.
#     }
#
#     if(disty%in%c(4,5,6)){
#       control$optiDisp <- as.integer(1)
#       theta.score <- as.numeric(rep(NA, length(theta)))
#     }else{
#       control$optiDisp <- as.integer(0)
#       theta.score <- -99999
#     }
#     scores <- as.numeric(rep(NA,length(c(alpha.score,beta.score,eta.score,
#                                          gamma.score,delta.score,theta.score))))
#     control$conv <- as.integer(0)
#
#     #model quantities
#     pis_out <- as.numeric(rep(NA, G))  #container for the fitted RCP model
#     mus <- as.numeric(array( NA, dim=c(n, S, G)))  #container for the fitted spp model
#     loglikeS <- as.numeric(rep(NA, S))
#     loglikeSG  <- as.numeric(matrix(NA, nrow = S, ncol = G))
#
#     optimiseDisp <- FALSE
#     kount <- 1
#     tmp.new <- tmp.old <- -999999
#     if( control$optimise.cpp){
#       while( (abs( abs( tmp.new - tmp.old) / ( abs( tmp.old) + control$reltol)) > control$reltol | kount==1) & (kount < 15)){
#         kount <- kount + 1
#         tmp.old <- tmp.new
#         message( "Updating Location Parameters: ", appendLF=FALSE)
#         tmp <- .Call("species_mix_cpp",
#                       as.numeric(as.matrix(y)), as.numeric(as.matrix(X)),
#                       as.numeric(as.matrix(Wcpp)), as.numeric(as.matrix(Ucpp)),
#                       as.numeric(offset), as.numeric(spp_weights),
#                       as.numeric(as.matrix(site_spp_weights)),
#                       as.integer(as.matrix(!y_is_na)),
#                       as.numeric(size), as.integer(S), as.integer(G), as.integer(npx),
#                       as.integer(npw), as.integer(npu), as.integer(n),
#                       as.integer(disty),as.integer(TRUE),
#                       as.integer(control$optiPart),as.integer(control$optiAll),
#                       # SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty, //data
#                       as.double(alpha), as.double(beta), as.double(eta),
#                       as.double(gamma), as.double(delta), as.double(theta),
#                       as.double(powers),
#                       # SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
#                       as.numeric(control$penalty.alpha),as.numeric(control$penalty.beta),
#                       as.numeric(control$penalty.pi),as.numeric(control$penalty.gamma),
#                       as.numeric(control$penalty.delta),
#                       as.numeric(control$penalty.theta[1]),
#                       as.numeric(control$penalty.theta[2]),
#                       # SEXP &RalphaPen, SEXP &RbetaPen, SEXP &RpiPen,  SEXP &RgammaPen,
#                       # SEXP &RdeltaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen,
#                       alpha.score, beta.score, eta.score, gamma.score, delta.score,
#                       theta.score, as.integer(control$getscores.cpp), as.numeric(scores),
#                       # SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
#                       pis_out, mus, loglikeS, loglikeSG,
#                       # SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
#                       as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport),
#                       as.numeric(control$abstol), as.numeric(control$reltol), as.integer(control$conv),
#                       as.integer(control$printparams.cpp),
#                       # SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
#                       as.integer(control$optimise.cpp), as.integer(TRUE), as.integer(FALSE),
#                       # SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp
#                       PACKAGE = "ecomix")
#         # tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
#                      # as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
#                      # alpha, tau, beta, gamma, disp, power,
#                      # as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
#                      # alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
#                      # pi, mus, logCondDens, logls,
#                      # as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
#                      # as.integer(control$optimise.cpp), as.integer(TRUE), as.integer( FALSE), as.integer(optimiseDisp), as.integer( FALSE), PACKAGE = "ecomix")
#         message( "Updating Dispersion Parameters: ", appendLF=FALSE)
#         for( ii in 1:S){
#           tmp1 <- nlminb( theta[ii], Tw.phi.func, Tw.phi.func.grad, spp3=ii, control=list( trace=0))
#           theta[ii] <- tmp1$par
#           message( tmp1$objective, " ")
#         }
#         message( "")
#         tmp.new <- -tmp1$objective
#       }
#     }
#     # tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
#     #              as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
#     #              alpha, tau, beta, gamma, disp, power,
#     #              as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
#     #              alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
#     #              pi, mus, logCondDens, logls,
#     #              as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
#     #              as.integer(FALSE), as.integer( TRUE), as.integer(TRUE), as.integer(TRUE), as.integer( FALSE), PACKAGE = "ecomix")
#
#     ret <- list()
#
#     ret$pi <- matrix(pi, ncol = nRCP)
#     ret$mus <- array( mus, dim=c(n,S,nRCP))
#     ret$coefs <- list(alpha = alpha, tau = tau, beta = beta, gamma=gamma, disp=disp)
#     if( any( ret$coefs$gamma==-999999, na.rm=TRUE))
#       ret$coefs$gamma <- NULL
#     if( any( ret$coefs$disp==-999999, na.rm=TRUE))
#       ret$coefs$disp <- NULL
#     ret$names <- list( spp=colnames( outcomes), RCPs=paste( "RCP", 1:nRCP, sep=""), Xvars=colnames( X))
#     if( p.w>0)
#       ret$names$Wvars <- colnames( W)
#     else
#       ret$names$Wvars <- NA
#     ret$scores <- list(alpha = alpha.score, tau = tau.score, beta = beta.score, gamma = gamma.score, disp=disp.score)
#     if( any( ret$scores$gamma==-999999, na.rm=TRUE))
#       ret$scores$gamma <- NULL
#     if( any( ret$scores$disp==-999999, na.rm=TRUE))
#       ret$scores$disp <- NULL
#     ret$logCondDens <- matrix(logCondDens, ncol = nRCP)
#     if( control$optimise)
#       ret$conv <- conv
#     else
#       ret$conv <- "not optimised"
#     ret$S <- S; ret$nRCP <- nRCP; ret$p.x <- p.x; ret$p.w <- p.w; ret$n <- n
#     ret$start.vals <- inits
#     ret$logl <- tmp
#     ret$logl.sites <- logls  #for residuals
#
#     return( ret)
#   }

"transform.theta" <- function(theta, disty, logify = TRUE){#, max.theta = NULL){

  if(logify){
  if(disty%in%4){
      log.new.theta <- log(theta)
  }
  if(disty%in%5){
    log.new.theta  <- log(theta)
  }
  if(disty%in%6){
    log.new.theta <- log(sqrt(theta))
  }
    return(log.new.theta)
  }  else {
    if(disty%in%4){
      new.theta <- exp(theta)
    }
    if(disty%in%5){
      new.theta  <- exp(theta)
    }
    if(disty%in%6){
      new.theta <- exp(theta)^2
    }
      return(new.theta)
    }
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

"calc_post_probs_sam" <- function( pi, logCondDens)  {
  logPostProbs <- log( pi) + logCondDens
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

"clean_data_sam" <- function(data, form1, form2, form3=NULL){#, family){
  # if(family=='ippm') na_rule <- "na.pass"
  # else
  na_rule <- "na.exclude"
  # form1 <- update.formula(form1, ~ . -1)
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

"get_family_sam" <- function(family, size) {

 if(class(family)=="family"){
    fam = family
    if(fam$family=="binomial" & fam$link=="cloglog")
      family="cloglog"
    if(fam$family=="binomial" & fam$link=="logit")
      family="binomial"
    if(fam$family=="poisson" & fam$link=="log")
      family="poisson"
    if(fam$family=="gaussian" & fam$link=="identity")
      family="gaussian"
  }

  if ( is.character(family) ) {
    if(family=="ippm")stop("'ippm' models are now run from the ecomix.ppm package\n")
    if (substr(family,1,3) == "ber") {
      disty <- 1  #logistic/binomial
      familyname <- "bernoulli"
      link <- "logit"  #by default
    }
    else if (substr(family,1,3) == "bin" & all(size==1)) {
      disty <- 1  #logistic/binomial
      familyname <- "bernoulli"
      link <- "logit"  #by default
    }
    else if (substr(family,1,3) == "bin" & !all(size==1)) {
      disty <- 7  #logistic/binomial
      familyname <- "binomial"
      link <- "logit"
    }
    else if (substr(family,1,3) == "clo") {
      disty <- 1  #binomial(cloglog) ## I might make this 3
      familyname <- "bernoulli"
      link <- "cloglog"
    }
    else if (substr(family,1,3) == "poi") {
      disty <- 2   #poisson
      familyname <- "poisson"
      link <- "log"
    }
    else if (substr(family,1,3) == "neg") {
      disty <- 4  #neg bin
      familyname <- "negative.binomial"
      link <- "log"
    }
    else if (substr(family,1,3) == "twe") {
      disty <- 5
      familyname<- "tweedie"
      link <- "log" # not meaningful, we only use the lo link at the moment
    }
    else if (substr(family,1,3) == "gau") {
      disty <- 6 # gaussian
      familyname <- "gaussian"
      link <- "idenity"
    }
    # FAMILY EDIT
    else stop (paste("'family'", family, "not recognised. See for currently available options."))
  }
  return(list(disty=disty,family=familyname,link=link))
}

"get_logls_sam" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
                            offset, y_is_na, disty, linky, size, powers, control,
                            fits, get_fitted=TRUE){

   if(get_fitted) fitted_values <- array(0,dim=c(G,nrow(y),S))

   #setup the right link function
   # if(is.null(linky)){
   # if(disty%in%c(1,7)) link <- stats::make.link(link = "logit")
   # if(disty%in%c(2,3,4,5)) link <- stats::make.link(link = "log")
   # if(disty%in%c(6)) link <- stats::make.link(link = "identity")
   # } else {
   link <- stats::make.link(link=linky)
   # }

   logl_sp <- matrix(NA, nrow=S, ncol=G)

   if(!is.null(U))eta_all <- as.matrix(U) %*% fits$delta
   else eta_all <- matrix(0,nrow=nrow(X),ncol=1)

    for(ss in 1:S){
      sp_idx<-!y_is_na[,ss]
      for(gg in 1:G){
        eta_mix <- X[sp_idx,,drop=FALSE] %*% fits$beta[gg,]
        if(ncol(W)>1){
          eta_spp <- W[sp_idx,,drop=FALSE] %*% c(fits$alpha[ss],fits$gamma[ss,])
        } else {
          eta_spp <- fits$alpha[ss]
        }
        eta <- eta_spp + eta_mix + eta_all[sp_idx,,drop=FALSE] + offset[sp_idx]
        if(get_fitted) fitted_values[gg,sp_idx,ss] <- link$linkinv(eta)

        if(disty==1) logl_sp[ss,gg] <- sum(dbinom(y[,ss], size =  1, prob =  link$linkinv(eta),log = TRUE))
        if(disty==2) logl_sp[ss,gg] <- sum(dpois(y[,ss], lambda = link$linkinv(eta),log = TRUE))
        # if(disty==3) logl_sp[ss,gg] <- y[sp_idx,ss] %*% eta - site_spp_weights[sp_idx,ss] %*% link$linkinv(eta)
        if(disty==4) logl_sp[ss,gg] <- sum(dnbinom(y[,ss], mu=link$linkinv(eta), size = 1/fits$theta[ss], log=TRUE))
        if(disty==5) logl_sp[ss,gg] <- sum(fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
        if(disty==6) logl_sp[ss,gg] <- sum(dnorm(y[,ss],mean=eta,sd=sqrt(fits$theta[ss]),log=TRUE))
        if(disty==7) logl_sp[ss,gg] <- sum(dbinom(y[,ss],size =  size, prob = link$linkinv(eta),log = TRUE))

      }
      # if(!disty%in%3)
      logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
    }

   out.list <- list(logl_sp=logl_sp)
   if(get_fitted) out.list$fitted = fitted_values
   return(out.list)
}

# "get_logls_spp_sam" <- function(y, X, W, U, G, S, spp_weights, site_spp_weights,
#                             offset, y_is_na, disty, linky, size, powers, control,
#                             fits){
#
#   #setup the right link function
#   # if(disty%in%c(1,7)) link <- stats::make.link(link = "logit")
#   # if(disty%in%c(2,3,4,5)) link <- stats::make.link(link = "log")
#   # if(disty%in%c(6)) link <- stats::make.link(link = "identity")
#   link <- stats::make.link(link=linky)
#
#   logl_sp <- matrix(NA, nrow=S, ncol=G)
#
#   if(!is.null(U))eta_all <- as.matrix(U) %*% fits$delta
#   else eta_all <- rep(0,nrow(X))
#
#   for(ss in 1:S){
#     sp_idx<-!y_is_na[,ss]
#     for(gg in 1:G){
#       eta_mix <- as.matrix(X[sp_idx,]) %*% fits$beta[gg,]
#       if(ncol(W)>1) eta_spp <- as.matrix(W[sp_idx,,drop=FALSE]) %*% c(fits$alpha[ss],fits$gamma[ss,])
#       else eta_spp <- fits$alpha[ss]
#       eta <- eta_spp + eta_mix + eta_all[sp_idx] + offset[sp_idx]
#
#       if(disty==1) logl_sp[ss,gg] <- sum(dbinom(y[,ss], 1, link$linkinv(eta),log = TRUE))
#       if(disty==2) logl_sp[ss,gg] <- sum(dpois(y[,ss], lambda = link$linkinv(eta),log = TRUE))
#       # if(disty==3) logl_sp[ss,gg] <- y[sp_idx,ss] %*% eta - site_spp_weights[sp_idx,ss] %*% link$linkinv(eta)
#       if(disty==4) logl_sp[ss,gg] <- sum(dnbinom(y[,ss], mu=link$linkinv(eta), size = fits$theta[ss], log=TRUE))
#       if(disty==5) logl_sp[ss,gg] <- sum(fishMod::dTweedie(y[,ss], mu =  link$linkinv(eta), phi = fits$theta[ss], p = powers[ss], LOG = TRUE))
#       if(disty==6) logl_sp[ss,gg] <- sum(dnorm(y[,ss],mean=eta,sd=sqrt(fits$theta[ss]),log=TRUE))
#       if(disty==7) logl_sp[ss,gg] <- sum(dbinom(y[,ss], size, link$linkinv(eta),log = TRUE))
#     }
#     # if(!disty%in%3)
#     logl_sp[ss,gg] <- logl_sp[ss,gg]*spp_weights[ss]
#   }
#
#   ak <- logl_sp + matrix(rep(log(fits$pi), each=S), nrow=S, ncol=G)
#   am <- apply( ak, 1, max)
#   ak <- exp( ak-am)
#   sppLogls <- am + log( rowSums( ak))
#   return(sppLogls)
#
# }

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


"get_site_spp_weights_sam"  <- function(mf, site_spp_weights, sp_names){#, family){
  if(is.null(site_spp_weights))site_spp_weights <- model.weights(mf)
  # if(family=='ippm'){
  #   if(!is.null(site_spp_weights)){
  #     site_spp_weights <- subset(site_spp_weights, select = colnames(site_spp_weights)%in%sp_names)
  #   } else {
  #     site_spp_weights <- matrix(1,nrow(mf),length(sp_names))
  #   }
  # } else {
    if(!is.null(site_spp_weights)){
      site_spp_weights <- replicate(length(sp_names),site_spp_weights)
    } else {
      site_spp_weights <- matrix(1,nrow(mf),length(sp_names))
    }
  # }
  return(site_spp_weights)
}

"get_titbits_sam" <- function(titbits, y, X, W, U, spp_weights, site_spp_weights, offset,
                              y_is_na, size, powers, archetype_formula, species_formula,
                              all_formula, data, control, family, link)  {
  if( titbits==TRUE)
    titbits <- list(Y = y, X = X, W = W, U = U, spp_weights = spp_weights,
                    site_spp_weights = site_spp_weights, offset = offset,
                    y_is_na = y_is_na, size = size, powers=powers,
                    archetype_formula =  archetype_formula,
                    species_formula = species_formula,
                    all_formula = all_formula, data = data, control = control,
                    family = family,link=link)
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
    if( "data" %in% titbits)
      titbits$data <- data
    if( "control" %in% titbits)
      titbits$control <- control
    if( "family" %in% titbits)
      titbits$family <- family
    if( "link" %in% titbits)
      titbits$link <- link
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

"delete.intercept" <- function(mm) {
  ## Save the attributes prior to removing the intercept coloumn:
  saveattr <- attributes(mm)
  ## Find the intercept coloumn:
  intercept <- which(saveattr$assign == 0)
  ## Return if there was no intercept coloumn:
  if (!length(intercept)) return(mm)
  ## Remove the intercept coloumn:
  mm <- mm[,-intercept, drop=FALSE]
  ## Update the attributes with the new dimensions:
  saveattr$dim <- dim(mm)
  saveattr$dimnames <- dimnames(mm)
  ## Remove the assignment of the intercept from the attributes:
  saveattr$assign <- saveattr$assign[-intercept]
  ## Restore the (modified) attributes:
  attributes(mm) <- saveattr
  ## Return the model matrix:
  mm
}


"get_X_sam" <- function(mf.X){

  mt.x <- terms(mf.X)
  mt.x <- stats::delete.response(mt.x)
  X <- stats::model.matrix(mt.x, mf.X)
  X <- delete.intercept(X)

  return(list(X=X,mt.x=mt.x))
}

"get_W_sam" <- function(mf.W){

  mt.w <- stats::terms(mf.W)
  mt.w <- stats::delete.response(mt.w)
  W <- stats::model.matrix(mt.w, mf.W)

  return(list(W=W, mt.w=mt.w))
}

"get_U_sam" <- function(mf.U){

  if(!is.null(mf.U)){
    mt.u <- terms(mf.U)
    mt.u <- stats::delete.response(mt.u)
    U <- stats::model.matrix(mt.u, mf.U)
  } else {
    U <- NULL
    mt.u <- NULL
  }
  return(list(U=U,mt.u=mt.u))
}

"print_input_sam" <- function(y, X, W, U, S, archetype_formula, species_formula,
                              all_formula, family, link, quiet=FALSE){
  if( quiet)
    return( NULL)
  n.tot <- nrow(y)
  message("There are ", nrow(X), " site observations for ", S," species")

  archetype_formula[[2]] <- NULL
  message("The model for the archetype (grouping) is ", Reduce( "paste", deparse(archetype_formula)))
  if(!is.null(species_formula))
    message("The model for the species is ", Reduce( "paste", deparse(species_formula)))
  if(!is.null(U))
    message("The model for the entire dataset is ", Reduce( "paste", deparse(all_formula)))
    message("You are implementing a ", family, " Species Archetype Model.")
    message("This model uses a ",link, " link function.")
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
  message("starting species specific dispersion parameters (these are on a log scale):\n",paste(round(theta,3)," "))
}

"sam_model_clean_up_names" <- function(mod){

  names(mod$coefs$alpha) <- mod$names$spp
  colnames(mod$coefs$beta) <- mod$names$Xvars
  rownames(mod$coefs$beta) <- mod$names$SAMs
  names(mod$pi) <- mod$names$SAMs
  if( any(mod$coefs$gamma==-999999, na.rm=TRUE)){
    mod$coefs$gamma <- NULL
  } else {
    colnames(mod$coefs$gamma) <- mod$names$Wvars
    rownames(mod$coefs$gamma) <- mod$names$spp
  }
  if( any(mod$coefs$delta==-999999, na.rm=TRUE)){
    mod$coefs$delta <- NULL
  } else {
    names(mod$coefs$delta) <- mod$names$Uvars
  }
  if( any( mod$coefs$theta==-999999, na.rm=TRUE)|any(!mod$disty%in%c(4,5,6))){
    mod$coefs$theta <- NULL
  } else {
    names(mod$coefs$theta) <- mod$names$spp
  }

  colnames(mod$tau) <- mod$names$SAMs
  rownames(mod$tau) <-mod$names$spp
  return(mod)

}


"sam_random_inits" <- function(alpha, beta, gamma, delta, theta, S, G, X, W, U, disty, mult=0.3, control.sd = control$init.sd){
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
  if(disty %in% c(4,5,6)){
    my.sd <- mult*sd( theta); if( is.na( my.sd) | my.sd==0) my.sd <- 0.1
    theta <- theta + as.numeric( rnorm( S, mean=0, my.sd))
  }
  if(disty %in% c(4,5,6)) return(list(alpha,beta,gamma,delta,theta))
  else return(list(alpha,beta,gamma,delta))
}


"sam_internal_pred_groups" <- function(alpha, beta, tau, gamma, delta,
                                       G, S, X, W, U, offset = NULL, family, type){

  if (family %in% c("bernoulli","binomial"))
    link.fun <- make.link("logit")
  if (family %in% c("negative.binomial","poisson","tweedie"))
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
      if(type=='response')s.outpred[, s] <- link.fun$linkinv(eta)
      else if (type=='link')s.outpred[, s] <- link.fun$linkinv(eta)
      else stop ('type not known')
    }

    outpred_arch[, g] <- apply(s.outpred*rep(tau[, g],each = dim(X)[1]),
                               1, sum)/sum(tau[, g])
  }
  return(outpred_arch)
}

"sam_internal_pred_species" <- function(alpha, beta, tau, gamma, delta,
                                        G, S, X, W, U, offset = NULL, family, type){

  if (family %in% c("bernoulli","binomial"))
    link.fun <- make.link("logit")
  if (family %in% c("negative.binomial","poisson","tweedie"))
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
    if(type=='response') mug <- link.fun$linkinv(eta)
    else if(type=='link') mug <- eta
    else stop('type not known')
    outpred_spp <- outpred_spp + mug*matrix(tau[,g], nrow(X), S, byrow=TRUE)
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
      theta <- rep(-999999,S)
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
  if (!("init.method" %in% names(control)))
    control$init.method <- 'kmeans'
  if (!("init.sd" %in% names(control)))
    control$init.sd <- NA
  if (!("minimum.sites" %in% names(control)))
    control$minimum.sites <- 0
  if (!("ecm.prefit" %in% names(control)))
    control$ecm.prefit <- TRUE # if(family!="ippm")
    # else control$ecm.prefit <- FALSE
  if (!("ecm.steps" %in% names(control)))
    control$ecm.steps <- 10
  if (!("ecm.refits" %in% names(control)))
    control$ecm.refits <- 1
  if (!("ecm.reltol" %in% names(control)))
    control$ecm.reltol <- 1e-3
  if (!("print.cpp.vals" %in% names(control)))
    control$print.cpp.vals <- FALSE
  if (!("nreport" %in% names(control)))
    control$nreport <- 10
  if (!("abstol" %in% names(control)))
    control$abstol <- 1e-05
  if (!("reltol" %in% names(control)))
    control$reltol <- sqrt(.Machine$double.eps)
  if (!("optimise.cpp" %in% names( control)))
    control$optimise.cpp <- TRUE
  if (!("getscores.cpp" %in% names( control)))
  control$getscores.cpp <- FALSE
  if (!("loglOnly.cpp" %in% names(control)))
    control$loglOnly.cpp <- FALSE
  if (!("derivOnly.cpp" %in% names( control)))
    control$derivOnly.cpp <- FALSE
  if (!("printparams.cpp" %in% names( control)))
    control$printparams.cpp <- FALSE
  if(!("doPenalties")%in%names(control))
    control$doPenalties <- 0
  if (!("penalty.pi" %in% names(control)))
    control$penalty.pi <- 0.01
  else
    if (control$penalty.pi < 0) {
      message("Supplied penalty for pi is negative, reverting to the default")
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

"shrink_taus" <- function(tau, G){
  if( G==1)
    return( tau)
  magical.alpha <- (1-0.8*G)/(0.8*(2-G)-1) ## Dunstan et al., 2013, JABES
  taus_star <- (2*magical.alpha*tau-magical.alpha+1)/(2*magical.alpha - magical.alpha*G + G)
  return(taus_star)
}

"species_data_check" <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))
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

