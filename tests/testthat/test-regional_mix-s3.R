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
