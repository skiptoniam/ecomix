test_that("species_mix.multifit with a single nArchetypes fits nstart models and prints a summary", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix.multifit(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                              data = d$data, family = "bernoulli", nArchetypes = 2, nstart = 2,
                              control = quiet_control)
  expect_s3_class(fm, "species_mix.multifit")
  expect_length(fm$multiple_fits, 2)
  expect_true(all(vapply(fm$multiple_fits, function(x) x$G == 2, logical(1))))
  expect_output(print(fm))
})

test_that("species_mix.multifit with group selection fits the correct nArchetypes at each position", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix.multifit(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                              data = d$data, family = "bernoulli", nArchetypes = c(2, 3), nstart = 2,
                              control = quiet_control)
  expect_s3_class(fm, "species_mix.multifit")
  expect_length(fm$multiple_fits, 2)
  expect_true(all(vapply(fm$multiple_fits[[1]], function(x) x$G == 2, logical(1))))
  expect_true(all(vapply(fm$multiple_fits[[2]], function(x) x$G == 3, logical(1))))
  expect_output(print(fm))
})

test_that("plot.species_mix.multifit doesn't turn a fully-degenerate nArchetypes column into Inf", {
  # Regression test: when every nstart refit for a given nArchetypes is
  # flagged as "ObviouslyBad" (an archetype collapsed to ~0 expected
  # species), min(x, na.rm=TRUE) on the resulting all-NA column used to
  # silently return Inf (with a warning) instead of NA, corrupting the
  # BIC/AIC/logLik curve -- see the SAM vignette's group-selection example.
  d <- make_sam_data("bernoulli")
  fm <- species_mix.multifit(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                              data = d$data, family = "bernoulli", nArchetypes = c(2, 3), nstart = 2,
                              control = quiet_control)
  # force every refit at the first nArchetypes value to look degenerate
  for (ii in seq_along(fm$multiple_fits[[1]])) {
    fm$multiple_fits[[1]][[ii]]$tau[, 1] <- 0
  }

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_warning(plot(fm, type = "BIC"))
})
