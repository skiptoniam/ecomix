test_that("effects_data.regional_mix and plot.regional_mix_effects_data work for continuous covariates", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  ed <- effects_data("x1", fm)
  expect_s3_class(ed, "regional_mix_effects_data")
  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    plot.regional_mix_effects_data(ed, fm)
  })
})
