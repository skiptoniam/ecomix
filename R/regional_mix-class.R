##### Main species mix functions to export #####
#' @title This is how you fit a regions of common profiles model in ecomix
#' @rdname regional_mix
#' @name regional_mix
#' @description This is how you fit \code{regional_mix} model in ecomix.
#' @param rcp_formula an object of class "formula" (or an object that can be coerced to that class). The response variable (left hand side of the formula) needs to be either 'presence', 'occurrence', 'abundance', 'biomass' or 'quantity' this will help specify the type of data to be modelled, if the response variable is disperate to the model family an error will be thrown. The dependent variables (the right hind side) of this formula specifies the dependence of the region of common profile (rcp) probabilities on covariates.
#' @param species_formula an object of class "formula" (or an object that can be coerced to that class). The left hand side of this formula should be left empty (it is removed if it is not empty). The right hand side of this formula specifies the dependence of the species"'" data on covariates (typically different covariates to \code{rcp_formula} to avoid confusing confounding). An example formula is observations ~ gear_type + time_of_day, where gear_type describes the different sampling gears and time_of_day describes the time of the sample. #maybe could call this detection/bias
#' @param data a List which contains named objects 'species_data': a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned. And 'covariate_data' a data frame containing the covariate data for each site. Names of columns must match that given in \code{rcp_formula} and \code{species_formula}.
#' @param nRCP The number of mixing components (groups) to fit.
#' @param family The family of statistical family to use within the ecomix models. a  choice between "bernoulli", "poisson", "negative.binomial", "tweedie" and "gaussian" distributions are possible and applicable to specific types of data.
#' @param offset a numeric vector of length nrow( data) that is included into the model as an offset. It is included into the conditional part of the model where conditioning is performed on the unobserved RCP type. Note that offsets cannot be included as part of the rcp_formula or species_formula arguments ??? only through this argument.
#' @param weights a numeric vector of length nrow( data) that is used as weights in the log-likelihood calculations. If NULL (default) then all weights are assumed to be identically 1.
#' @param control a list of control parameters for optimisation and calculation. See details. From \code{control} control.
#' @param inits a character string which defines the method used to initialise finite mixture model clustering. #Will have to synergise this function call across RCP and SpeciesMix. Looks like SpeciesMix uses a em.prefit to setup initialisations. regional_mix has a number of methods. This seems like a good place to setup the bivariate clustering step - cobra function.
#' @param titbits either a boolean or a vector of characters. If TRUE (default for regional_mix(qv)), then some objects used in the estimation of the model"'"s parameters are returned in a list entitled "titbits" in the model object. Some functions, for example plot.regional_mix(qv) and predict.regional_mix(qv), will require some or all of these pieces of information. If titbits=FALSE (default for regional_mix.multifit(qv)), then an empty list is returned. If a character vector, then just those objects are returned. Possible values are:"Y" for the outcome matrix, "X" for the model matrix for the RCP model, "W" for the model matrix for the species-specific model, "offset" for the offset in the model, "wts" for the model weights, "form.RCP" for the formula for the RCPs, "form.spp" for the formula for the species-specific model, "control" for the control arguments used in model fitting, "family" for the conditional family of the species data, and "power" for the power parameters used (only used in Tweedie models). Care needs to be taken when using titbits=TRUE in regional_mix.multifit(qv) calls as titbits is created for EACH OF THE MODEL FITS. If the data is large or if nstart is large, then setting titbits=TRUE may give users problems with memory. It is more efficient, from a memory perspective, to refit the "best" model using regional_mix(qv) after identifying it with regional_mix.multifit(qv). See examples for illustration about how to do this.
#' @param power a numeric vector (length either 1 or the number of species) defining the power parameter to use in the Tweedie models. If length(power)==1, then the same power parameter is used for all species. If length(power)==No_species, then each species gets its own power parameter. Power values must be between 1 and 2, for computational reasons they should be well away from the boundary. The default is 1.6 as this has proved to be a good ball-park value for the fisheries data that the developer has previously analysed.
#' @importFrom graphics abline hist legend lines matplot par plot points polygon rect
#' @importFrom stats as.formula binomial cooks.distance cov cutree dbinom dist dnbinom dnorm dpois
#' fitted gaussian glm hclust lm logLik model.matrix model.offset model.response
#' model.weights pbinom pnbinom pnorm poisson
#' ppois predict qnorm qqnorm quantile rbinom
#' residuals rgamma rnbinom rnorm rpois runif
#' sd uniroot update update.formula
#' @details A typical formula for use in the rcp_formula argument will have the
#' form (for example) cbind(spp1,spp2,spp3,spp4)~1+cov1+cov2*cov3. This signifies
#' that there are 4 species to be used for RCP modelling and that the RCP types
#' are dependent on cov1+cov2+cov3+cov2:cov3. See ?glm for a description of how
#' the right hand side of the formula is expanded.
#'
#' Likewise a typical formula for use in the species_formula argument will have
#' the form (for example) ~1+fac1+cov1. This signifies that the catchabilities
#' of each species depends upon the levels of the factor fac1 and the covariate
#' cov1. See ?glm for a description of how the right hand side of the formula
#' is expanded.
#'
#' The computation strategy for the default method, which has been demonstrated
#' to work for all data sets the developers have encountered thus far, is fully
#' described in Foster et al (2013) and Foster et al (2017). We note however,
#' that it is a good idea to standardise covariates prior to calling regional_mix.
#' This is not formally required by the model, but it does drastically reduce
#' the chance of numerical issues in the first iteration. If you choose to NOT
#' standardise, then you should at least choose a scale that is reasonable
#' (so that the numerical range is measured by units and not thousands of units).
#' This may mean that the units may be, for example, kilometres (and not metres),
#' or 100s of kilometres (and not metres/kilometres).
#'
#' We do not, on purpose, provide residuals as a routine part of the model.
#' Users should use the residuals.regional_mix(qv) function to obtain them. We do
#' this as the type of residual needs to be specified (although we recommend
#' type=="RQR" for routine use).
#'
#' Control arguments for optimisation generally follow those in optim(qv),
#' although a few differences occur (e.g. "loglOnly").
#' The elements of the control list are
#' \describe{
#'  \item{maxit}{The maximum number of iterations. Default is 500.}
#'  \item{quiet}{Should any reporting be performed? Default is FALSE, for reporting. For regional_mix.multifit(), this indicates if the progress should be printed.}
#'  \item{trace}{Non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information.}
#'  \item{nreport}{The frequency of reports for optimisation. Default is 10 ??? a report for 10th iteration.}
#'  \item{reltol}{Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to sqrt(.Machine$double.eps), typically about 1e-8.}
#'  \item{optimise}{Should optimisation for estimation occur? If TRUE (default) optimisation will occur. If FALSE no optimisation is performed.}
#'  \item{loglOnly}{Should the log-likelihood be caulcated? If TRUE (default) then log-likelihood is calculated and returned. If FALSE then the log-likelihood is not calculated for return.}
#'  \item{derivOnly}{Should the scores be evaluated at the (final) parameter values. If TRUE (default) then they are calculated. If FALSE then they are not calculated.}
#'  \item{penalty}{A numeric scalar. This is the concentration for the Dirichlet-inspired penalty for the prior probabilities. Values less than zero will be set to the default (0.1). Large values give more penalisation than small ones.}
#'  \item{penalty.tau}{A numeric scalar. This is the penalty for the tau parameters in the species model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
#'  \item{penalty.gamma}{A numeric scalar. This is the penalty for the gamma parameters in the species model. They are assumed to come from a normal distribution with standard deviation given as this parameter (default is 10).}
#'  \item{penalty.disp}{a two element vector. These are combined to form the penalty for the dispersion parameters (if any). The dispersions are assumed to come from a log-normal distribution with log-mean penalty.disp[1] and log-standard-deviation penalty.disp[2]. Defaults to c(10,sqrt(10)), which gives shrinkage towards 1 (the mode of the penalty). Note that for Normal models, where the dispersion alone defines the variance, a strong penalty may be required to keep parameters estimable.}
#'  }
#'  For calls to regimix.multifit(), titbits is set to FALSE??? so no excess memory is used. If users want this information, and there is good reason to want it, then a call to regimix() with starting values given as the best fit's estimates should be used.
#' @return regional_mix returns an object of class \code{regional_mix} and regional_mix.multifit returns a list of objects of class \code{regional_mix}. The regional_mix class has several methods: coef, plot, predict, residuals, summary, and vcov. The regional_mix object consists of a list with the following elements:
#' @return AIC Akaike an information criterion for the maximised model.
#' @return BIC Bayesian information criterion for the maximised model.
#' @return call the call to the function.
#' @return coefs a list of three elements, one each for the estimates for the species prevalence (alpha), the deviations from alpha for the first (nRCP-1) RCP (tau), and the (nRCP-1) sets of RCP regression coefficients (beta).
#' @return conv the convergence code from the maximisation procedure. See ?optim for an explanation (basically 0 is good and anything else is bad).
#' @return dist the character string identifying the family used for the model.
#' @return logCondDens an nObs by nRCP matrix specifying the probability of observing each sites"'" data, given each of the RCP types.
#' @return logl the maximised log likelihood.
#' @return mus an array of size nRCP x S x nRCP where each element of the first dimension is the fitted value for all the species in all the RCP types.
#' @return n the number of samples.
#' @return names the names of the species, and the names of the covariates for X and W.
#' @return nRCP the number of RCPs.
#' @return pis an n x nRCP matrix with each column giving the prior probabilities for the corresponding RCP type. Rows sum to one.
#' @return postProbs an n x nRCP matrix with each column giving the posterior probabilities for the corresponding RCP type. Rows sum to one (as each site is assumed to be from one of the RCP types).
#' @return p.w the number of covariates used in the species-specific model.
#' @return p.x the number of covariates used in the RCP model
#' @return S the number of species.
#' @return scores a list of three elements. Structure corresponds to coefs.
#' @return start.vals the values used to start the estimation procedure.
#' @return titbits (if requested using the titbit argument, see above) other pieces of information, useful to developers, that users should not typically need to concern themselves with. However, this information is used by methods for regional_mix objects.
#' @export
#' @examples
#' \dontrun{
#' rcp_form <- as.formula(paste0("cbind(",paste(colnames(simulated_data[,1:20]),
#' collapse = ','),")~1+x1+x2+x3"))
#' spp_form <- observations ~ 1 + w1 + w2
#' fm_regional_mix <- regional_mix(rcp_formula=rcp_form,species_formula=spp_form,
#'                                 data=data, family='bernoulli', nRCP=5)
#' }
"regional_mix" <- function (rcp_formula = NULL, species_formula = NULL, data,
                            nRCP = 3, family="bernoulli", offset=NULL,
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
  Xdat <- get_X_rcp(dat$mf.X) #rcp_formula,
  X <- Xdat$X
  xterms <- Xdat$mt.x
  p.x <- ncol( X)
  #get design matrix for spp part of the model -- if there is one
  Wdat <- get_W_rcp( species_formula, dat$mf.W)
  W <- Wdat$W
  wterms <- Wdat$mt.w
  if( all( W != -999999))
    p.w <- ncol( W)
  else
    p.w <- 0
  #get offset (if not specified then it will be zeros)
  offy <- get_offset_rcp( mf, dat$mf.X, dat$mf.W)
  #get model wts (if not specified then it will be ones)
  wts <- get_wts_rcp( mf)
  #get family
  disty.cases <- c("bernoulli","poisson","negative.binomial","tweedie",
                   "gaussian")
  disty <- get_dist_rcp( disty.cases, family)

  ## catch the terms as a list
  tt <- list(xterms=xterms,wterms=wterms)

  #get power params for Tweedie
  power <- get_power_rcp( disty, power, S)
  #summarising data to console
  print.data.summ( data, dat, S, rcp_formula, species_formula, disty.cases,
                   disty, control$quiet)

  tmp <- regional_mix.fit( outcomes, W, X, offy, wts, disty, nRCP, power, inits,
                           control, nrow( X), S, p.x, p.w)

  tmp$family <- disty.cases[disty]
  #calculate the posterior probs
  if( nRCP>1)
    tmp$postProbs <- calcPostProbs( tmp$pis, tmp$logCondDens)
  else
    tmp$postProbs <- rep( 1, nrow( X))
  #Residuals --not calculating residuals here.
  #Information criteria
  tmp <- calcInfoCrit( tmp)
  #titbits object, if wanted/needed.
  tmp$titbits <- get_titbits_rcp( titbits, outcomes, X, W, offy, wts, data,
                                  rcp_formula, species_formula, control,
                                  family, p.w=p.w, power)
  tmp$titbits$disty <- disty
  #the last bit of the regional_mix object puzzle
  tmp$call <- call
  tmp$terms <- tt
  tmp <- tmp[sort( names( tmp))]
  gc()
  class(tmp) <- "regional_mix"
  return(tmp)

  #documentation needs to be adjusted to fit new model.

}

#'@title regional_mix.fit
#'@rdname regional_mix.fit
#'@name regional_mix.fit
#' @description regional_mix.fit is similar to glm.fit and does all the heavy lifting when it
#' comes to estimating regional mix models. If you are unfamilar with how to use glm.fit it is
#' recommended that you use regional_mix which is the user friendly wrapper around this function.
#'@param outcomes is a matrix genertated from model.response containing the species information. The matrix has the dimensions n_sites * n_species.
#'@param W is a design matrix for regional_formula and will be implemented if regional_formula has covariates.
#'@param X is a design matrix for the archetype_formula dimension n_sites * n_covariates.
#'@param offy this is a vector of site specific offsets, this might be something like the log(area sampled at sites).
#'@param wts is the site weights. These are weights used to alter the loglikelihood.
#'@param disty the error family to used in regional_mix estimation. Currently, 'bernoulli', 'poisson', 'negative.binomial' and 'guassian' are available - internal conversion of family to a integer.
#'@param nRCP is the number of species archetypes that are being estimated.
#'@param control this is a list of control parameters that alter the specifics of model fitting.
#'@param power This is for the Tweedie distribution - currently not in use (until we fix the Tweedie computational stuff).
#'@param inits This will be a vector of starting values for regional_mix (i.e you've fitted a model and want to refit it).
#'@param n the number of sites in the data
#'@param S is the number of species to be modelled (this will be calculated internally in regional_mix())
#'@param p.x The number of covariates fitted to the design matrix X.
#'@param p.w The number of covariates fitted to the design matrix W.
#'@export

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
    if( disty != 4){ #not Tweedie
    optimiseDisp <- TRUE
    tmp <- notTweedieOptimise( outcomes, X, W, offy, wts, S, nRCP, p.x,
                               p.w, nrow( X), disty, start.vals, power, control)
    } else { #Tweedie -- quite convoluted in comparison
    tmp <- TweedieOptimise( outcomes, X, W, offy, wts, S, nRCP, p.x, p.w, nrow( X), disty, start.vals, power, control)
    }
    return( tmp)

  }

#' @rdname regional_mix
#' @name regional_mix.multifit
#' @param nstart for regional_mix.multifit only. The number of random starts to perform for re-fitting. Default is 10, which will need increasing for serious use.
#' @param mc.cores for regional_mix.multifit only. The number of cores to spread the re-fitting over.
#' @export

"regional_mix.multifit" <-  function(rcp_formula = NULL, species_formula = NULL,
                                     data, nRCP = 3, family="bernoulli",
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
    Xdat <- get_X_rcp(dat$mf.X) #rcp_formula,
    X <- Xdat$X
    xterms <- Xdat$mt.x
    p.x <- ncol( X)
    #get design matrix for spp part of the model -- if there is one
    Wdat <- get_W_rcp( species_formula, dat$mf.W)
    W <- Wdat$W
    wterms <- Wdat$mt.w
    if( all( W != -999999))
      p.w <- ncol( W)
    else
      p.w <- 0

    ## catch the terms as a list
    tt <- list(xterms=xterms,wterms=wterms)

    #get offset (if not specified then it will be zeros)
    offy <- get_offset_rcp( mf, dat$mf.X, dat$mf.W)
    #get model wts (if not specified then it will be ones)
    wts <- get_wts_rcp( mf)
    #get family
    disty.cases <- c("bernoulli","poisson",
                     "negative.binomial","tweedie","gaussian")
    disty <- get_dist_rcp( disty.cases, family)
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
      tmp$family <- disty.cases[disty]
      #calculate the posterior probs
      if( nRCP>1)
        tmp$postProbs <- calcPostProbs( tmp$pis, tmp$logCondDens)
      else
        tmp$postProbs <- rep( 1, nrow( X))
      #Residuals --not calculating residuals here.  Need to call residuals.regional_mix
      #Information criteria
      tmp <- calcInfoCrit( tmp)
      #titbits object, if wanted/needed.
      tmp$titbits <- get_titbits_rcp( titbits, outcomes, X, W, offy, wts, data,
                                      rcp_formula, species_formula, control,
                                      family, p.w=p.w, power)
      tmp$titbits$disty <- disty
      #the last bit of the regional_mix object puzzle
      tmp$call <- call
      tmp$terms <- tt
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

#### Non S3 Class objects ####
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

"check_if_sampling" <-function(object){object$p.w<1}

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
  error.msg <- paste( c( "family not implemented. Options are: ", disty.cases, "-- Exitting Now"), collapse=" ")
  disty <- switch( dist1, "bernoulli" = 1,"poisson" = 2,"negative.binomial" = 3,"tweedie" = 4,"gaussian" = 5,{stop( error.msg)} )
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
function( site.logls, outcomes, family, coef, nRCP, type="deviance", powers=NULL, quiet=FALSE, nsim=1000, X, W, offy)
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
    sims <- regional_mix.simulate( nRCP=nRCP, S=length( coef$alpha), n=n.sim, p.x=ncol( X), p.w=ncol( W), alpha=coef$alpha, tau=coef$tau, beta=coef$beta, gamma=coef$gamma, logDisps=coef$disp, powers=pwers, X=X1, W=W1, offset=offy,family=family)

  }

  return( resids)

}


"get_start_vals_rcp" <- function( outcomes, W, X, offy, wts, disty, G, S, power, inits, quiet=FALSE){
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
function( titbits, outcomes, X, W, offset, wts, data, rcp_formula, species_formula, control, family, p.w, power)
{
  if( titbits==TRUE)
    titbits <- list( Y = outcomes, X = X, W = W, offset = offset, wts=wts, data =data, rcp_formula = rcp_formula, species_formula = species_formula, control = control, family = family, power=power)
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
    if( "data" %in% titbits)
      titbits$data <- data
    if( "rcp_formula" %in% titbits)
      titbits$rcp_formula <- rcp_formula
    if( "species_formula" %in% titbits)
      titbits$species_formula <- species_formula
    if( "control" %in% titbits)
      titbits$control <- control
    if( "family" %in% titbits)
      titbits$family <- family
    if( "power" %in% titbits)
      titbits$power <- power
  }
  if( p.w==0 & "W" %in% names( titbits))
    titbits$W <- NULL
  if( p.w!=0 & "species_formula" %in% names( titbits))
    environment( titbits$species_formula) <- environment( titbits$rcp_formula)
  return( titbits)

}


# "get_W_rcp" <- function( species_formula, mf.W){
#   form.W <- species_formula
#   if( !is.null( species_formula)){
#     if( length( form.W)>2)
#       form.W[[2]] <- NULL #get rid of outcomes
#     W <- model.matrix( form.W, mf.W)
#     tmp.fun <- function(x){ all( x==1)}
#     intercepts <- apply( W, 2, tmp.fun)
#     W <- W[,!intercepts,drop=FALSE]
#   }
#   else
#     W <- -999999
#   return( W)
# }

# "get_W_rcp" <- function (form.spp, mf.W){
#   form.W <- form.spp
#   if (!is.null(form.spp)) {
#     if (length(form.spp) > 2)
#       form.W[[2]] <- NULL
#     W <- model.matrix(form.W, mf.W)
#     tmp.fun <- function(x) {
#       all(x == 1)
#     }
#     intercepts <- apply(W, 2, tmp.fun)
#     W <- W[, !intercepts, drop = FALSE]
#   }
#   else W <- -999999
#   return(W)
# }

"get_W_rcp" <- function(form.spp,mf.W){

  if (!is.null(form.spp)) {
    # if (length(form.spp) > 2)
      # form.W[[2]] <- NULL
    mt.w <- stats::terms(mf.W)
    mt.w <- stats::delete.response(mt.w)
    mt.w <- delete.intercept(mt.w)
    W <- stats::model.matrix(mt.w, mf.W)
  } else {
    W <- -999999
    mt.w <- -999999
  }
  return(list(W=W, mt.w=mt.w))
}



"get_wts_rcp" <-function ( mf){
  wts <- model.weights( mf)
  if( is.null( wts))
    return( rep( 1, nrow( mf)))  #all weights assumed equal
  return( wts)
}

"get_X_rcp" <-function( mf.X){ #rcp_formula,
  # old code.
  # form.X <- rcp_formula
  # form.X[[2]] <- NULL
  # form.X <- as.formula(form.X)
  # X <- model.matrix(form.X, mf.X)

  mt.x <- terms(mf.X)
  mt.x <- stats::delete.response(mt.x)
  mt.x <- delete.response(mt.x)
  X <- stats::model.matrix(mt.x, mf.X)

  return( list(X=X,mt.x=mt.x))
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
    if( disty==3){  #negative.binomial
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
  message("The error family is: ", disty.cases[disty])
}

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


"TweedieOptimise" <-
function( outcomes, X, W, offy, wts, S, nRCP, p.x, p.w, n, disty, start.vals, power, control){
  Tw.phi.func <- function( phi1, spp3){
    disp3 <- disp
    disp3[spp3] <- phi1
    tmp1 <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
      as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
      alpha, tau, beta, gamma, disp3, power,
      as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
      alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
      pis, mus, logCondDens, logls,
      as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
      as.integer(FALSE), as.integer(TRUE), as.integer( FALSE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")
    return( -as.numeric( tmp1))
  }

 Tw.phi.func.grad <- function( phi1, spp3){
   disp3 <- disp
   disp3[spp3] <- phi1
   tmp.disp.score <- rep( -99999, S)
   tmp1 <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
     as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
     alpha, tau, beta, gamma, disp3, power,
     as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
     alpha.score, tau.score, beta.score, gamma.score, tmp.disp.score, scoreContri,
     pis, mus, logCondDens, logls,
     as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
     as.integer(FALSE), as.integer(FALSE), as.integer(TRUE), as.integer( TRUE), as.integer( FALSE), PACKAGE = "ecomix")
   return( -as.numeric( tmp.disp.score[spp3]))
 }

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
 scoreContri <- -999999  #as.numeric(matrix(NA, ncol = length(inits), nrow = n))
 #model quantities
 pis <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))  #container for the fitted RCP model
 mus <- as.numeric(array( NA, dim=c( n, S, nRCP)))  #container for the fitted spp model
 logCondDens <- as.numeric(matrix(NA, nrow = n, ncol = nRCP))
 logls <- as.numeric(rep(NA, n))
 conv <- as.integer(0)

 optimiseDisp <- FALSE
 kount <- 1
 tmp.new <- tmp.old <- -999999
 if( control$optimise){
   while( (abs( abs( tmp.new - tmp.old) / ( abs( tmp.old) + control$reltol)) > control$reltol | kount==1) & (kount < 15)){
     kount <- kount + 1
     tmp.old <- tmp.new
     message( "Updating Location Parameters: ", appendLF=FALSE)
     tmp <- .Call("RCP_C", as.numeric(outcomes), as.numeric(X), as.numeric(W), as.numeric( offy), as.numeric( wts),
       as.integer(S), as.integer(nRCP), as.integer(p.x), as.integer(p.w), as.integer(n), as.integer( disty),
       alpha, tau, beta, gamma, disp, power,
       as.numeric(control$penalty), as.numeric(control$penalty.tau), as.numeric( control$penalty.gamma), as.numeric( control$penalty.disp[1]), as.numeric( control$penalty.disp[2]),
       alpha.score, tau.score, beta.score, gamma.score, disp.score, scoreContri,
       pis, mus, logCondDens, logls,
       as.integer(control$maxit), as.integer(control$trace), as.integer(control$nreport), as.numeric(control$abstol), as.numeric(control$reltol), as.integer(conv),
       as.integer(control$optimise), as.integer(TRUE), as.integer( FALSE), as.integer(optimiseDisp), as.integer( FALSE), PACKAGE = "ecomix")
     message( "Updating Dispersion Parameters: ", appendLF=FALSE)
     for( ii in 1:S){
       tmp1 <- nlminb( disp[ii], Tw.phi.func, Tw.phi.func.grad, spp3=ii, control=list( trace=0))
       disp[ii] <- tmp1$par
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



