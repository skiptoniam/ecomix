% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional_mix-class.R
\name{regional_mix.fit}
\alias{regional_mix.fit}
\title{regional_mix.fit}
\usage{
regional_mix.fit(
  outcomes,
  W,
  X,
  offy,
  wts,
  disty,
  nRCP,
  power,
  inits,
  control,
  n,
  S,
  p.x,
  p.w
)
}
\arguments{
\item{outcomes}{is a matrix genertated from model.response containing the species information. The matrix has the dimensions n_sites * n_species.}

\item{W}{is a design matrix for regional_formula and will be implemented if regional_formula has covariates.}

\item{X}{is a design matrix for the archetype_formula dimension n_sites * n_covariates.}

\item{offy}{this is a vector of site specific offsets, this might be something like the log(area sampled at sites).}

\item{wts}{is the site weights. These are weights used to alter the loglikelihood.}

\item{disty}{the error family to used in regional_mix estimation. Currently, 'bernoulli', 'poisson', 'negative.binomial' and 'guassian' are available - internal conversion of family to a integer.}

\item{nRCP}{is the number of species archetypes that are being estimated.}

\item{power}{This is for the Tweedie distribution - currently not in use (until we fix the Tweedie computational stuff).}

\item{inits}{This will be a vector of starting values for regional_mix (i.e you've fitted a model and want to refit it).}

\item{control}{this is a list of control parameters that alter the specifics of model fitting.}

\item{n}{the number of sites in the data}

\item{S}{is the number of species to be modelled (this will be calculated internally in regional_mix())}

\item{p.x}{The number of covariates fitted to the design matrix X.}

\item{p.w}{The number of covariates fitted to the design matrix W.}
}
\description{
regional_mix.fit is similar to glm.fit and does all the heavy lifting when it
comes to estimating regional mix models. If you are unfamilar with how to use glm.fit it is
recommended that you use regional_mix which is the user friendly wrapper around this function.
}
