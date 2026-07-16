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
