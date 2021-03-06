% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regional_mix_s3-class.R,
%   R/species_mix_s3-class.R
\name{effectPlotData.regional_mix}
\alias{effectPlotData.regional_mix}
\alias{effectPlotData}
\alias{effectPlotData.species_mix}
\title{Generate data for plotting or predicting partial effects of covariates}
\usage{
\method{effectPlotData}{regional_mix}(focal.predictors, mod, ngrid = 50, ...)

effectPlotData(focal.predictors, mod, ...)

\method{effectPlotData}{species_mix}(focal.predictors, mod, ngrid = 50, ...)
}
\arguments{
\item{focal.predictors}{A character or string of characters which represent covariates in the model.}

\item{mod}{The fitted species_mix model.}

\item{ngrid}{The length of the prediction vector.}

\item{...}{other arguments}
}
\value{
This function should return a list of data.frames one for each focal.predictor.
This will enable user to predict marginal effects or plot the partial response plots.
}
\description{
This function produces a list of data.frames for predicting the partial
effect of a focal.predictor current included in a species_mix model.
}
\details{
Generate data for plotting or predicting partial effects of covariates

This function produces a list of data.frames for predicting the partial
effect of a focal.predictor current included in a species_mix model.
}
\examples{
\donttest{
library(ecomix)
set.seed(42)
sam_form <- stats::as.formula(paste0('cbind(',paste(paste0('spp',1:20),
collapse = ','),")~x1+x2"))
sp_form <- ~ 1
beta <- matrix(c(-2.9,-3.6,-0.9,1,.9,1.9),3,2,byrow=TRUE)
dat <- data.frame(y=rep(1,100),x1=stats::runif(100,0,2.5),
x2=stats::rnorm(100,0,2.5))
dat[,-1] <- scale(dat[,-1])
simulated_data <- species_mix.simulate(archetype_formula = sam_form,species_formula = sp_form,
data = dat,beta=beta,family="bernoulli")
fm1 <- species_mix(archetype_formula = sam_form,species_formula = sp_form,
data = simulated_data, family = 'bernoulli',  nArchetypes=3)
effectPlotData("x1",fm1)
}
}
