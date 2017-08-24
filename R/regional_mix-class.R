#' @useDynLib ecomix
# #' @importFrom Rcpp sourceCpp
NULL

#' @title ecomix package
#' @name ecomix-package
#' @description ecomix is a package for implementing species archetype models (SAM) and region of
#' common profile models (RCP) in R. The main function is \code{ecomix}, which fits either a SAM or RCP
#' based on input parameters.
#' @docType package
NULL

#' @title regional_mix objects
#' @rdname regional_mix
#' @name regional_mix
#' @description creates an \code{regional_mix} model.
#' @param rcp_formula an object of class "formula" (or an object that can be coerced to that class).
#' The response variable (left hand side of the formula) needs to be either 'presence', 'occurrence', 'abundance', 'biomass' or 'quantity' this will help specify the type of data to be modelled, if the response variable is disperate to the model distribution an error will be thrown. The dependent variables (the right hind side) of this formula specifies the dependence of the region of common profile (rcp) probabilities on covariates.
#' @param species_formula an object of class "formula" (or an object that can be coerced to that class). The left hand side of this formula should be left empty (it is removed if it is not empty). The right hand side of this formula specifies the dependence of the species"'" data on covariates (typically different covariates to \code{rcp_formula} to avoid confusing confounding). An example formula is observations ~ gear_type + time_of_day, where gear_type describes the different sampling gears and time_of_day describes the time of the sample. #maybe could call this detection/bias
#' @param model_data a List which contains named objects 'species_data': a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned. And 'covariate_data' a data frame containng the covariate data for each site. Names of columns must match that given in \code{rcp_formula} and \code{species_formula}.
#' @param n_mixtures The number of mixing components (groups) to fit.
#' @param distribution The family of statistical distribution to use within the ecomix models. a  choice between "bernoulli", "poisson", "negative_binomial", "tweedie" and "gaussian" distributions are possible and applicable to specific types of data.
#' @param offset a numeric vector of length nrow( data) that is included into the model as an offset. It is included into the conditional part of the model where conditioning is performed on the unobserved RCP type. Note that offsets cannot be included as part of the rcp_formula or species_formula arguments ??? only through this argument.
#' @param weights a numeric vector of length nrow( data) that is used as weights in the log-likelihood calculations. If NULL (default) then all weights are assumed to be identically 1.
#' @param control a list of control parameters for optimisation and calculation. See details. From \code{control} control.
#' @param initialise a characture string which defines the method used to initialise finite mixture model clustering. #Will have to synergise this function call across RCP and SpeciesMix. Looks like SpeciesMix uses a em.prefit to setup initialisations. regional_mix has a number of methods. This seems like a good place to setup the bivariate clusting step - cobra function.
#' @export
#' @examples
#' simulated_data <- simulate_regional_mix_data()
#' rcp_form <- as.formula(paste0("cbind(",paste(sort(sp_name),collapse = ','),")~1+x1+x2+x3"))
#' spp_form <- observations ~ 1 + w1 + w2
#' model_data <- list('species_data' = Y, 'covariate_data' = X)
#' fm_regional_mix <- regional_mix(rcp_form,spp_form,data=model_data,distribution='bernoulli',n_mixtures=5)

regional_mix <- function(rcp_formula=NULL,
                       species_formula=NULL,
                       model_data,
                       n_mixtures=3,
                       distribution = 'bernoulli',
                       initialise='random',
                       offset=NULL,
                       weights=NULL,
                       control=list(),
                       titbits=TRUE,
                       power=1.6){
  stop('haha this does nothing yet.')

  #use this to match formula to distribution in model
  #presence -> poisson w/ offset
  #occurrence -> bernoulli
  #abundance -> poisson or negative_binomial
  #biomass -> tweedie
  #quantity -> guassian
  #all.vars(process_formula)[1]

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
    m <- match(c("model_data","offset","weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.exclude"
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ##model_data <- as.data.frame(model_data)
    #get the data model frames and strip out any NAs
    dat <- clean.data( mf, rcp_formula, species_formula)
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
    X <- get.X(rcp_formula, dat$mf.X)
    p.x <- ncol( X)
    #get design matrix for spp part of the model -- if there is one
    W <- get.W( species_formula, dat$mf.W)
    if( all( W != -999999))
      p.w <- ncol( W)
    else
      p.w <- 0
    #get offset (if not specified then it will be zeros)
    offy <- get.offset( mf, dat$mf.X, dat$mf.W)
    #get model wts (if not specified then it will be ones)
    wts <- get.wts( mf)
    #get distributionribution
    disty.cases <- c("Bernoulli","Poisson","NegBin","Tweedie","Normal")
    disty <- get.dist( disty.cases, distribution)
    #get power params for Tweedie
    power <- get.power( disty, power, S)
    #summarising data to console
    print.data.summ( model_data, dat, S, rcp_formula, species_formula, disty.cases, disty, control$quiet)

    tmp <- regimix.fit( outcomes, W, X, offy, wts, disty, nRCP, power, inits, control, nrow( X), S, p.x, p.w)

    tmp$distribution <- disty.cases[disty]
    #calculate the posterior probs
    if( nRCP>1)
      tmp$postProbs <- calcPostProbs( tmp$pis, tmp$logCondDens)
    else
      tmp$postProbs <- rep( 1, nrow( X))
    #Residuals --not calculating residuals here.  Need to call residuals.regimix
    #Information criteria
    tmp <- calcInfoCrit( tmp)
    #titbits object, if wanted/needed.
    tmp$titbits <- get.titbits( titbits, outcomes, X, W, offy, wts, rcp_formula, species_formula, control, distribution, p.w=p.w, power)
    tmp$titbits$disty <- disty
    #the last bit of the regimix object puzzle
    tmp$call <- call

    gc()
    tmp <- tmp[sort( names( tmp))]

    class(tmp) <- "regimix"
    return(tmp)

    #documentation needs to be adjusted to fit new model.

  }
