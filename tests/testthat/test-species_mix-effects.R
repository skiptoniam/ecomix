test_that("effects_data.species_mix and plot.species_mix_effects_data work for continuous covariates", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control,
                     titbits = TRUE)
  ed <- effects_data("x1", fm)
  expect_s3_class(ed, "species_mix_effects_data")
  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    plot(ed, fm)
  })
})

test_that("effects_data.species_mix varies the focal predictor and holds others at their mean", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control,
                     titbits = TRUE)
  ed <- effects_data(focal.predictors = c("x1", "x2"), fm, ngrid = 10)
  expect_named(ed, c("x1", "x2"))
  expect_equal(range(ed$x1$x1), range(fm$titbits$data$x1))
  expect_length(unique(ed$x1$x2), 1)
  expect_equal(unique(ed$x1$x2), mean(fm$titbits$data$x2), ignore_attr = TRUE)
  expect_length(unique(ed$x2$x1), 1)
})

test_that("plot.species_mix_effects_data reports real bootstrap uncertainty when given a boot.object", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control,
                     titbits = TRUE)
  set.seed(2)
  sam_boot <- bootstrap(fm, nboot = 6, quiet = TRUE)
  ed <- effects_data("x1", fm)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm, boot.object = sam_boot))
})

test_that("plot.species_mix_effects_data works for Species, SpeciesSum and a single named species", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control,
                     titbits = TRUE)
  set.seed(2)
  sam_boot <- bootstrap(fm, nboot = 6, quiet = TRUE)
  ed <- effects_data("x1", fm)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm, boot.object = sam_boot, response.var = "Species"))
  expect_no_error(plot(ed, fm, boot.object = sam_boot, response.var = "SpeciesSum"))
  # regression test: subsetting object$tau to a single species used to drop
  # to a vector and crash apply(..., 1, which.max)
  expect_no_error(plot(ed, fm, boot.object = sam_boot, response.var = fm$names$spp[1]))
  expect_no_error(plot(ed, fm, response.var = fm$names$spp[1]))
})

test_that("plot.species_mix_effects_data works for a partial SAM (species_formula with covariates)", {
  d <- make_sam_data_partial("bernoulli", seed = 42)
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control,
                     titbits = TRUE)
  ed <- effects_data(c("x1", "x3"), fm)
  expect_named(ed, c("x1", "x3"))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm))
})
