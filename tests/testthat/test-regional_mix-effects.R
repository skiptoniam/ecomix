test_that("effects_data.regional_mix and plot.regional_mix_effects_data work for continuous covariates", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  ed <- effects_data("x1", fm)
  expect_s3_class(ed, "regional_mix_effects_data")
  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    plot(ed, fm)
  })
})

test_that("effects_data.regional_mix varies the focal predictor over its observed range", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  ed <- effects_data("x1", fm, ngrid = 10)
  expect_equal(range(ed$x1$x1), range(fm$titbits$X[, "x1"]))
  expect_length(ed$x1$x1, 10)
})

test_that("plot.regional_mix_effects_data reports real bootstrap uncertainty when given object2", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  set.seed(2)
  rcp_boot <- bootstrap(fm, nboot = 6, quiet = TRUE)
  ed <- effects_data("x1", fm)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm, object2 = rcp_boot))
})

test_that("plot.regional_mix_effects_data works for a subset of RCPs and a single named RCP", {
  d <- make_rcp_data("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, data = d$data, family = "bernoulli",
                      nRCP = 2, control = quiet_control)
  set.seed(2)
  rcp_boot <- bootstrap(fm, nboot = 6, quiet = TRUE)
  ed <- effects_data("x1", fm)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm, object2 = rcp_boot, response.var = "RCP_1"))
  expect_no_error(plot(ed, fm, response.var = "RCP_1"))
})

test_that("plot.regional_mix_effects_data works with a species_formula (p.w > 0)", {
  d <- make_rcp_data_partial("bernoulli")
  fm <- regional_mix(rcp_formula = d$rcp_formula, species_formula = d$species_formula,
                      data = d$data, family = "bernoulli", nRCP = 2, control = quiet_control)
  ed <- effects_data(c("x1", "w1"), fm)
  expect_named(ed, c("x1", "w1"))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm))
})

test_that("plot.regional_mix_effects_data works for a factor covariate, with and without a bootstrap object", {
  set.seed(1)
  n <- 80; S <- 6; nRCP <- 2
  x1 <- runif(n, -2, 2)
  grp <- factor(sample(c("a", "b"), n, replace = TRUE))
  Xmat <- model.matrix(~x1 + grp)
  colnames(Xmat) <- c("(Intercept)", "x1", "grpb")
  form.RCP <- as.formula(paste0("cbind(", paste0("spp", 1:S, collapse = ","), ")~x1+grp"))
  sim <- suppressMessages(regional_mix.simulate(nRCP = nRCP, S = S, p.x = 3, p.w = 0, n = n,
    alpha = rnorm(S, 1, 0.3),
    tau = matrix(rnorm((nRCP - 1) * S, 0, 0.5), nrow = nRCP - 1),
    beta = matrix(rnorm((nRCP - 1) * 3, 0, 0.3), nrow = nRCP - 1),
    X = Xmat, family = "bernoulli"))
  dat <- cbind(sim, x1 = x1, grp = grp)
  fm <- regional_mix(rcp_formula = form.RCP, data = dat, family = "bernoulli", nRCP = 2,
                      control = quiet_control)

  ed <- effects_data("grp", fm)
  expect_true(is.factor(ed$grp$grp))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(ed, fm))

  set.seed(2)
  rcp_boot <- bootstrap(fm, nboot = 6, quiet = TRUE)
  expect_no_error(plot(ed, fm, object2 = rcp_boot))
})
