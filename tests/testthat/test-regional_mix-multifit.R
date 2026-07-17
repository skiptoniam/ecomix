test_that("regional_mix.multifit fits nstart models with class and a working print method", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix.multifit(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                               nRCP = 2, nstart = 2, control = quiet_control)
  expect_s3_class(fm, "regional_mix.multifit")
  expect_length(fm$multiple_fits, 2)
  expect_output(print(fm))
})

test_that("check_RCP_posteriors flags an artificially tiny RCP", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  expect_warning(res <- check_RCP_posteriors(fm, min_membership = 1000))
  expect_equal(sum(res), fm$n)

  fm.multi <- regional_mix.multifit(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                                     nRCP = 2, nstart = 2, control = quiet_control)
  res.multi <- suppressWarnings(check_RCP_posteriors(fm.multi, min_membership = 1000))
  expect_length(res.multi, 2)
})
