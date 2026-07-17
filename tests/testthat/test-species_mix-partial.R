test_that("partial species_mix (species_formula with covariates) fits and moves gamma off its starting values", {
  d <- make_sam_data_partial("bernoulli", seed = 42)
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.finite(fm$logl))
  expect_true(!is.null(fm$coefs$gamma))
  expect_equal(dim(fm$coefs$gamma), c(8, 1))

  npw_actual <- ncol(fm$coefs$gamma)
  start_idx <- fm$S + (fm$G * fm$npx) + (fm$G - 1) + seq_len(fm$S * npw_actual)
  start_gamma <- fm$start.vals[start_idx]
  expect_false(isTRUE(all.equal(as.numeric(fm$coefs$gamma), as.numeric(start_gamma))))

  expect_gt(stats::cor(as.numeric(fm$coefs$gamma), as.numeric(d$gamma)), 0.9)
})

