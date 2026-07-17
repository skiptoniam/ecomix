test_that("AIC/BIC/logLik/extractAIC return finite scalars consistent with the fitted logl", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  expect_true(is.finite(AIC(fm)))
  expect_true(is.finite(BIC(fm)))
  expect_equal(unclass(logLik(fm)), fm$logl)
  expect_gt(AIC(fm), -2 * fm$logl)
  expect_gt(BIC(fm), -2 * fm$logl)

  eaic <- extractAIC(fm)
  expect_length(eaic, 2)
  expect_equal(eaic[1], length(unlist(coef(fm))))
  expect_equal(eaic[2], AIC(fm))
})

test_that("coef.regional_mix returns correctly shaped alpha/tau/beta", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  co <- coef(fm)
  expect_named(co$alpha, fm$names$spp)
  expect_equal(dim(co$tau), c(fm$nRCP - 1, fm$S))
  expect_equal(dim(co$beta), c(fm$nRCP - 1, fm$p.x))
  expect_null(co$gamma)
})

test_that("coef.regional_mix returns dispersion for negative.binomial models", {
  d <- make_rcp_data("negative.binomial")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "negative.binomial",
                      nRCP = 2, control = quiet_control)
  co <- coef(fm)
  expect_named(co$logDisp, fm$names$spp)
})

test_that("print.regional_mix runs without error and returns its coefficient list invisibly", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  out <- expect_output(ret <- print(fm), "Distribution")
  expect_named(ret, c("Call", "Distribution", "coef"))
})

test_that("residuals.regional_mix returns finite RQR residuals of the right shape", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  resids <- residuals(fm, quiet = TRUE)
  expect_equal(dim(resids), c(fm$n, fm$S))
  expect_true(mean(is.finite(resids)) > 0.9)
})

test_that("residuals.regional_mix (deviance) returns one residual per site", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  resids <- suppressMessages(residuals(fm, type = "deviance", quiet = TRUE))
  expect_length(resids, fm$n)
})

test_that("vcov.regional_mix works with FiniteDifference and BayesBoot", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)

  vfd <- vcov(fm, method = "FiniteDifference")
  expect_true(is.matrix(vfd))
  expect_equal(dim(vfd), rep(length(unlist(fm$coefs)), 2))
  expect_true(isTRUE(all.equal(vfd, t(vfd))))
  expect_true(all(is.finite(vfd)))

  set.seed(1)
  vbb <- vcov(fm, method = "BayesBoot", nboot = 8)
  expect_true(is.matrix(vbb))
  expect_equal(dim(vbb), dim(vfd))
  expect_true(all(is.finite(vbb)))
})

test_that("summary.regional_mix returns an Estimate/SE/z-score/p table once vcov is attached", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  fm$vcov <- vcov(fm, method = "FiniteDifference")
  s <- expect_message(summary(fm), "Standard errors")
  expect_equal(colnames(s), c("Estimate", "SE", "z-score", "p"))
  expect_equal(nrow(s), length(unlist(fm$coefs)))
  expect_true(all(s[, "SE"] > 0))
})

test_that("summary.regional_mix errors informatively without a vcov", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  expect_error(summary(fm), "No variance matrix")
})

test_that("plot.regional_mix (RQR) runs without error", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(fm, type = "RQR"))
})

test_that("regional_mix.species_profile returns an RCP x species matrix without a species_formula", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  prof <- regional_mix.species_profile(fm)
  expect_s3_class(prof, "regional_mix_profile")
  expect_equal(dim(prof), c(fm$nRCP, fm$S))
  expect_true(all(prof >= 0 & prof <= 1))
})

test_that("regional_mix.species_profile gives a distinct profile per species-formula level", {
  d <- make_rcp_data_partial("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, species_formula = d$species_formula,
                      data = d$data, family = "bernoulli", nRCP = 2, control = quiet_control)
  prof <- regional_mix.species_profile(fm)
  expect_named(prof, c("(Intercept)", "w1"))
  expect_equal(dim(prof[["(Intercept)"]]), c(fm$nRCP, fm$S))
  # regression test: the gamma (species-covariate) effect must actually be
  # applied, not silently dropped -- levels should differ
  expect_false(isTRUE(all.equal(prof[["(Intercept)"]], prof[["w1"]])))
})

test_that("regional_mix.species_profile with a bootstrap object gives finite, varying summaries", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  set.seed(1)
  rcp_boot <- bootstrap(fm, nboot = 8, quiet = TRUE)
  prof <- regional_mix.species_profile(fm, rcp_boot)
  expect_named(prof, c("mean", "sd", "lower", "upper"))
  expect_equal(dim(prof$mean), c(fm$nRCP, fm$S))
  expect_true(all(is.finite(unlist(prof))))
  expect_true(all(prof$sd > 0))
  expect_true(all(prof$lower <= prof$mean & prof$mean <= prof$upper))
})

test_that("regional_mix.species_profile with a bootstrap object and species_formula gives distinct, varying levels", {
  d <- make_rcp_data_partial("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, species_formula = d$species_formula,
                      data = d$data, family = "bernoulli", nRCP = 2, control = quiet_control)
  set.seed(1)
  rcp_boot <- bootstrap(fm, nboot = 8, quiet = TRUE)
  prof <- regional_mix.species_profile(fm, rcp_boot)
  expect_named(prof, c("(Intercept)", "w1", "overall"))
  expect_true(all(is.finite(unlist(prof))))
  expect_true(all(prof[["(Intercept)"]]$sd > 0))
  expect_true(all(prof[["w1"]]$sd > 0))
  expect_false(isTRUE(all.equal(prof[["(Intercept)"]]$mean, prof[["w1"]]$mean)))
})

test_that("plot.regional_mix (deviance) runs without error", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control, titbits = TRUE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressMessages(plot(fm, type = "deviance", nsim = 5, quiet = TRUE)))
})

test_that("vcov.regional_mix works with SimpleBoot and EmpiricalInfo", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)

  set.seed(1)
  vsb <- vcov(fm, method = "SimpleBoot", nboot = 8)
  expect_true(is.matrix(vsb))
  expect_equal(dim(vsb), rep(length(unlist(fm$coefs)), 2))
  expect_true(all(is.finite(vsb)))

  vei <- suppressMessages(vcov(fm, method = "EmpiricalInfo"))
  expect_true(is.matrix(vei))
  expect_equal(dim(vei), dim(vsb))
  expect_true(all(is.finite(vei)))
})

test_that("residuals.regional_mix (RQR) works for poisson, negative.binomial and gaussian", {
  for (fam in c("poisson", "negative.binomial", "gaussian")) {
    d <- make_rcp_data(fam)
    fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = fam,
                        nRCP = 2, control = quiet_control)
    resids <- residuals(fm, quiet = TRUE)
    expect_equal(dim(resids), c(fm$n, fm$S), info = fam)
    expect_true(mean(is.finite(resids)) > 0.9, info = fam)
  }
})

test_that("regional_mix.simulate generates random parameters when none are supplied, and handles an unknown family", {
  expect_message(sim <- regional_mix.simulate(nRCP = 2, S = 5, n = 30, p.x = 2, p.w = 0,
                                               family = "bernoulli"),
                  "Random")
  expect_equal(dim(sim), c(30, 5 + 2))
  expect_true(all(as.matrix(sim[, 1:5]) %in% c(0, 1)))

  expect_message(res <- regional_mix.simulate(nRCP = 2, S = 5, n = 20, p.x = 2, p.w = 0,
                                               family = "notafamily"),
                  "family not found")
  expect_true(is.na(res))
})
