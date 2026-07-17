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

test_that("predict.species_mix with a bootstrap object gives species-level uncertainty (sam_cpp_pred species path)", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  set.seed(2)
  boot.obj <- bootstrap(fm, nboot = 8, type = "BayesBoot", mc.cores = 1, quiet = TRUE)

  p <- predict(fm, boot.object = boot.obj, prediction.type = "species", mc.cores = 1)
  expect_named(p, c("ptPreds", "bootPreds", "bootSEs", "bootCIs"))
  expect_equal(dim(p$ptPreds), c(fm$n, fm$S))
  expect_equal(dim(p$bootPreds), dim(p$ptPreds))
  expect_equal(dim(p$bootSEs), dim(p$ptPreds))
  expect_equal(dim(p$bootCIs), c(dim(p$ptPreds), 2))
  expect_true(all(p$ptPreds >= 0 & p$ptPreds <= 1))
  expect_true(all(p$bootSEs >= 0))
  expect_true(any(p$bootSEs > 0))
  expect_true(all(p$bootCIs[, , 1] <= p$bootCIs[, , 2]))

  # species-level bootstrap predictions also need to work with newdata, and
  # for a partial SAM (species_formula with real covariates)
  dp <- make_sam_data_partial("bernoulli", seed = 42)
  fmp <- species_mix(archetype_formula = dp$archetype_formula, species_formula = dp$species_formula,
                      data = dp$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  set.seed(3)
  boot.objp <- bootstrap(fmp, nboot = 8, type = "BayesBoot", mc.cores = 1, quiet = TRUE)
  newdata <- fmp$titbits$data[1:10, ]
  p2 <- predict(fmp, boot.object = boot.objp, newdata = newdata, prediction.type = "species", mc.cores = 1)
  expect_equal(dim(p2$ptPreds), c(10, fmp$S))
  expect_true(all(p2$bootSEs >= 0))
  expect_true(any(p2$bootSEs > 0))
})
