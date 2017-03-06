#' @useDynLib ecomix
#' @importFrom Rcpp sourceCpp
NULL

#' @title ecomix package
#' @name ecomix-package
#' @description ecomix is a package for implementing species archetype models (SAM) and region of
#' common profile models (RCP) in R. The main function is \code{ecomix}, which fits either a SAM or RCP
#' based on input parameters.
#' @docType package
NULL

#' @title ecomix objects
#' @rdname ecomix
#' @name ecomix
#' @description creates an \code{ecomix} model; see \code{ecomix.fit} for more specifics about how
#' ecomix fits finite mixture models.
#' @param process_formula an object of class "formula" (or an object that can be coerced to that class).
#' The response variable (left hand side of the formula) needs to be either 'presence', 'occurrence', 'abundance', 'biomass' or 'quantity' this will help specify the type of data to be modelled, if the response variable is disperate to the model distribution an error will be thrown. The dependent variables (the right hind side) of this formula specifies the dependence of the ecomix probabilities on covariates.
#' @param observation_formula an object of class "formula" (or an object that can be coerced to that class). The left hand side of this formula should be left empty (it is removed if it is not empty). The right hand side of this formula specifies the dependence of the species"'" data on covariates (typically different covariates to \code{process_formula} to avoid confusing confounding). An example formula is observations ~ gear_type + time_of_day, where gear_type describes the different sampling gears and time_of_day describes the time of the sample. #maybe could call this detection/bias
#' @param species_data a data frame containing the species information. The frame is arranged so that each row is a site and each column is a species. Species names should be included as column names otherwise numbers from 1:S are assigned.
#' @param covariate_data a data frame containng the covariate data for each site. Names of columns must match that given in \code{process_formula} and \code{observation_formula}.
#' @param n_mixtures The number of mixing components (groups) to fit.
#' @param distribution The family of statistical distribution to used within the ecomix models. a  choice between "bernoulli", "poisson", "negative_binomial", "tweedie" and "gaussian" distributions are possible and applicable to specific types of data.
#' @param offset a numeric vector of length nrow( data) that is included into the model as an offset. It is included into the conditional part of the model where conditioning is performed on the unobserved RCP type. Note that offsets cannot be included as part of the form.RCP or form.spp arguments ??? only through this argument.
#' @param weights a numeric vector of length nrow( data) that is used as weights in the log-likelihood calculations. If NULL (default) then all weights are assumed to be identically 1.
#' @param  control a list of control parameters for optimisation and calculation. See details. From \code{ecomix.fit} control.
#' @param initialise a characture string which defines the method used to initialise finite mixture model clustering. #Will have to synergise this function call across RCP and SpeciesMix. Looks like SpeciesMix uses a em.prefit to setup initialisations. regimix has a number of methods. This seems like a good place to setup the bivariate clusting step - cobra function.
#' @export
#' @examples
#' simulated_data <- simulate_ecomix_data()
#' process_form <- occurrence ~ 1 + x1 + x2 + x3
#' obs_form <- observations ~ 1 + w1 + w2 + w3

ecomix <- function(process_formula=NULL,
                   observation_formula=NULL,
                   species_data,
                   covariate_data,
                   n_mixtures=3,
                   initialise='random',
                   offset=NULL,
                   weights=NULL,
                   control=list()){
  message('haha this does nothing yet.')

  #use this to match formula to distribution in model
  #presence -> poisson w/ offset
  #occurrence -> bernoulli
  #abundance -> poisson or negative_binomial
  #biomass -> tweedie
  #quantity -> guassian
  #all.vars(process_formula)[1]



}
