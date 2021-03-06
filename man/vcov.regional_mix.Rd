% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional_mix_s3-class.R
\name{vcov.regional_mix}
\alias{vcov.regional_mix}
\title{Variance matrix for a regional_mix object.}
\usage{
\method{vcov}{regional_mix}(
  object,
  ...,
  object2 = NULL,
  method = "FiniteDifference",
  nboot = 1000,
  mc.cores = 1,
  D.accuracy = 2
)
}
\arguments{
\item{object}{an object obtained from fitting a RCP (for region of common profile) mixture model. Such as that generated from a call to regional_mix(qv).}

\item{object2}{an object of class \code{regional_mix} containing bootstrap samples of the parameter estimates (see regional_mix.bootstrap(qv)). If NULL (default) the bootstrapping is performed from within the vcov function. If not null, then the vcov estimate is obtained from these bootstrap samples.}

\item{method}{the method to calculate the variance-covariance matrix. Options are:'FiniteDifference' (default), \code{BayesBoot}, \code{SimpleBoot}, and \code{EmpiricalInfo}. The two bootstrap methods (\code{BayesBoot} and \code{SimpleBoot}, see regional_mix.bootstrap(qv)) should be more general and may possibly be more robust. The \code{EmpiricalInfo} method implements an empirical estimate of the Fisher information matrix, I can not recommend it however. It seems to behave poorly, even in well behaved simulations. It is computationally thrifty though.}

\item{nboot}{the number of bootstrap samples to take for the bootstrap estimation. Argument is ignored if !method \%in\% c(\code{FiniteDifference},'EmpiricalInfo').}

\item{mc.cores}{the number of cores to distrbute the calculations on. Default is 4. Set to 1 if the computer is running Windows (as it cannot handle forking -- see mclapply(qv)). Ignored if method=='EmpiricalInfo'.}

\item{D.accuracy}{The number of finite difference points used in the numerical approximation to the observed information matrix. Ignored if method != \code{FiniteDifference}. Options are 2 (default) and 4.}

\item{\\dots}{Other calls to the vcov function.}
}
\value{
A square matrix of size equal to the number of parameters. It contains the variance matrix of the parameter estimates.
}
\description{
Calculates variance-covariance matrix from a regional_mix object
}
\details{
If method is \code{FiniteDifference}, then the estimates variance matrix is based on a finite difference approximation to the observed information matrix.
If method is either "BayesBoot" or "SimpleBoot", then the estimated variance matrix is calculated from bootstrap samples of the parameter estimates. See Foster et al (in prep) for details of how the bootstrapping is actually done, and regional_mix.bootstrap(qv) for its implementation.
}
\references{
Foster, S.D., Givens, G.H., Dornan, G.J., Dunstan, P.K. and Darnell, R. (2013) Modelling Regions of Common Profiles Using Biological and Environmental Data. Environmetrics \emph{24}: 489--499. DOI: 10.1002/env.2245
Foster, S.D., Hill, N.A. and Lyons, M., 2017. Ecological grouping of survey sites when sampling artefacts are present. Journal of the Royal Statistical Society: Series C (Applied Statistics), 66(5), pp.1031-1047.
}
