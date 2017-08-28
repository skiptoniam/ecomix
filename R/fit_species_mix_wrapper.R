#this is a function which should take new forumula and data stuctures and fit to existing SAM data structure.
#old version of species mix took the formula obs ~ 1 + x1 + x2
#new version of species mix takes cbind(sp1, sp2, sp3) ~ 1 + x1 + x2
#so the trick will be setting up the model.frame and responses to suit the old approach.

fit_species_mix_wrapper <- function(formula, y, X, weights, offset, distribution, n_mixtures, inits, control, estimate_variance){

  sp.form <- update(form,obs~1+.)
  sp.data <- y
  covar.data <- X
  G <- n_mixtures
  pars <- init
  em.prefit <- control$em.prefit
  em.steps <- control$em.steps
  em.refit <- control$em.refit
  est.var <- control$estimate_variance
  trace <- control$trace
  r1 <- control$r1

  if(dist=="bernoulli") return(SpeciesMix.bernoulli(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace,r1))
  if(dist=="negative_binomial") return(SpeciesMix.nbinom(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace))
  if(dist=="tweedie") return(SpeciesMix.tweedie(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace))
  if(dist=="gaussian") return(SpeciesMix.gaussian(sp.form,sp.data,covar.data,G, pars, em.prefit,em.steps, em.refit ,est.var,residuals,trace))
  print("incorrect distribution type, options are bernoulli, negbin or tweedie")
}
