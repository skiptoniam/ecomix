test_that("AIC/BIC/logLik return finite scalars consistent with the fitted logl", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  expect_true(is.finite(AIC(fm)))
  expect_true(is.finite(BIC(fm)))
  expect_equal(unclass(logLik(fm)), fm$logl)
  # both are -2*logl plus a positive penalty on the number of free parameters
  expect_gt(AIC(fm), -2 * fm$logl)
  expect_gt(BIC(fm), -2 * fm$logl)
})

test_that("coef.species_mix returns correctly shaped alpha/beta for a non-partial model", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  co <- coef(fm)
  expect_named(co$alpha, fm$names$spp)
  expect_equal(dim(co$beta), c(fm$G, fm$npx))
  expect_null(co$gamma)
  expect_null(co$theta)
})

test_that("coef.species_mix returns a species x covariate gamma matrix for a partial SAM", {
  d <- make_sam_data_partial("bernoulli", seed = 42)
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  co <- coef(fm)
  expect_equal(dim(co$gamma), c(fm$S, fm$npw))
  expect_equal(fm$npw, ncol(d$gamma))
})

test_that("coef.species_mix returns dispersion for negative.binomial models", {
  d <- make_sam_data("negative.binomial")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "negative.binomial", nArchetypes = 2, control = quiet_control)
  co <- coef(fm)
  expect_named(co$theta, fm$names$spp)
  expect_true(all(co$theta > 0))
})

test_that("print.species_mix and terms.species_mix run without error", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  expect_output(print(fm), "species_mix model")
  tt <- terms(fm)
  expect_named(tt, c("xterms", "wterms"))
})

test_that("species_membership returns posterior taus that sum to 1 across archetypes", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  sm <- species_membership(fm)
  expect_s3_class(sm, "species_membership")
  expect_equal(dim(sm), c(fm$S, fm$G))
  expect_equal(unname(rowSums(sm)), rep(1, fm$S), tolerance = 1e-6)
})

test_that("residuals.species_mix returns finite RQR residuals of the right shape", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  resids <- residuals(fm, quiet = TRUE)
  expect_equal(dim(resids), c(fm$n, fm$S))
  expect_true(mean(is.finite(resids)) > 0.9)
})

test_that("vcov.species_mix works with FiniteDifference and BayesBoot for a non-partial model", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)

  vfd <- vcov(fm, method = "FiniteDifference")
  expect_true(is.matrix(vfd))
  expect_equal(dim(vfd), c(length(unlist(fm$coefs)), length(unlist(fm$coefs))))
  expect_true(isTRUE(all.equal(vfd, t(vfd))))
  expect_true(all(is.finite(vfd)))

  set.seed(1)
  vbb <- vcov(fm, method = "BayesBoot", nboot = 8)
  expect_true(is.matrix(vbb))
  expect_equal(dim(vbb), dim(vfd))
  expect_true(all(is.finite(vbb)))
})

test_that("vcov.species_mix (BayesBoot) runs for a partial SAM without error", {
  d <- make_sam_data_partial("bernoulli", seed = 1)
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  set.seed(1)
  v <- vcov(fm, method = "BayesBoot", nboot = 8)
  expect_true(is.matrix(v))
  expect_equal(dim(v), rep(length(unlist(fm$coefs)), 2))
})

test_that("summary.species_mix returns an Estimate/SE/z-score/p table once vcov is attached", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  fm$vcov <- vcov(fm, method = "FiniteDifference")
  s <- summary(fm)
  expect_equal(colnames(s), c("Estimate", "SE", "z-score", "p"))
  expect_equal(nrow(s), length(unlist(fm$coefs)))
  expect_true(all(s[, "SE"] > 0))
})

test_that("summary.species_mix errors informatively without a vcov", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  expect_error(summary(fm), "No variance matrix")
})
