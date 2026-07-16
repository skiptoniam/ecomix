test_that("regional_mix fits bernoulli data and returns a sane object", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  expect_s3_class(fm, "regional_mix")
  expect_true(is.finite(fm$logl))
  expect_equal(sum(fm$coefs$tau == 0), 0)  # tau (species-level RCP logit offsets) should be estimated
})

test_that("regional_mix fits poisson data", {
  d <- make_rcp_data("poisson")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "poisson",
                      nRCP = 2, control = quiet_control)
  expect_s3_class(fm, "regional_mix")
  expect_true(is.finite(fm$logl))
})

test_that("regional_mix fits negative.binomial data", {
  d <- make_rcp_data("negative.binomial")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "negative.binomial",
                      nRCP = 2, control = quiet_control)
  expect_s3_class(fm, "regional_mix")
  expect_true(is.finite(fm$logl))
  expect_true(!is.null(fm$coef$disp))
})

test_that("regional_mix fits gaussian data", {
  d <- make_rcp_data("gaussian")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "gaussian",
                      nRCP = 2, control = quiet_control)
  expect_s3_class(fm, "regional_mix")
  expect_true(is.finite(fm$logl))
})

test_that("regional_mix fits tweedie data", {
  d <- make_rcp_data("tweedie")
  fm <- suppressMessages(regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "tweedie",
                      power = 1.6, nRCP = 2, control = quiet_control))
  expect_s3_class(fm, "regional_mix")
  expect_true(is.finite(fm$logl))
})
