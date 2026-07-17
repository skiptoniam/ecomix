test_that("predict.species_mix gives archetype and species level point predictions", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  p_arch <- predict(fm, prediction.type = "archetype")
  p_spp <- predict(fm, prediction.type = "species")
  expect_equal(dim(p_arch), c(nrow(d$data), 2))
  expect_equal(dim(p_spp), c(nrow(d$data), 6))
  expect_true(all(p_arch >= 0 & p_arch <= 1))
  expect_true(all(p_spp >= 0 & p_spp <= 1))
})

test_that("bootstrap.species_mix produces varying dispersion estimates across replicates", {
  d <- make_sam_data("negative.binomial")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "negative.binomial", nArchetypes = 2, control = quiet_control)
  boot.obj <- bootstrap(fm, nboot = 4, type = "BayesBoot", mc.cores = 1, quiet = TRUE)
  expect_s3_class(boot.obj, "species_mix.bootstrap")
  expect_equal(nrow(boot.obj), 4)

  S <- fm$S; G <- fm$G; npx <- fm$npx; npw <- fm$npw
  thetaBoot <- boot.obj[, S + (G - 1) + (G * npx) + (S * npw) + seq_len(S), drop = FALSE]
  expect_equal(ncol(thetaBoot), S)
  # theta must vary across bootstrap replicates, not be frozen at a fixed value
  expect_true(all(apply(thetaBoot, 2, sd) > 0))
})

test_that("predict.species_mix with a bootstrap object gives bootPreds/bootSEs/bootCIs that reflect real uncertainty", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  boot.obj <- bootstrap(fm, nboot = 4, type = "BayesBoot", mc.cores = 1, quiet = TRUE)
  p <- predict(fm, boot.object = boot.obj, prediction.type = "archetype", mc.cores = 1)
  expect_named(p, c("ptPreds", "bootPreds", "bootSEs", "bootCIs"))
  expect_true(all(p$bootSEs >= 0))
  expect_true(any(p$bootSEs > 0))
})
