test_that("predict.regional_mix gives sane point predictions", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  p <- predict(fm)
  expect_equal(dim(p), c(nrow(d$data), 2))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(rowSums(p), rep(1, nrow(p)), tolerance = 1e-6)
})

test_that("predict.regional_mix with a bootstrap object gives bootPreds/bootSEs/bootCIs", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  boot.obj <- bootstrap(fm, nboot = 4, type = "BayesBoot", mc.cores = 1, quiet = TRUE)
  expect_s3_class(boot.obj, "regional_mix.bootstrap")
  p <- predict(fm, object2 = boot.obj, mc.cores = 1)
  expect_named(p, c("ptPreds", "bootPreds", "bootSEs", "bootCIs"))
  expect_true(all(p$bootSEs >= 0))
})
