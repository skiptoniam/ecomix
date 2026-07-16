make_sam_data <- function(family = "bernoulli", nArchetypes = 2, S = 6, n = 80,
                           seed = 1, power = 1.6) {
  set.seed(seed)
  sam_form <- stats::as.formula(paste0("cbind(",
    paste(paste0("spp", 1:S), collapse = ","), ")~x1+x2"))
  sp_form <- ~1
  beta <- matrix(rnorm(nArchetypes * 2, sd = 1), nrow = nArchetypes, ncol = 2)
  dat <- data.frame(y = rep(1, n), x1 = stats::runif(n, 0, 2.5),
                     x2 = stats::rnorm(n, 0, 2.5))
  dat[, -1] <- scale(dat[, -1])
  sim <- suppressMessages(species_mix.simulate(archetype_formula = sam_form,
    species_formula = sp_form, data = dat, beta = beta,
    nArchetypes = nArchetypes, powers = if (family == "tweedie") rep(power, S) else NULL,
    family = family))
  list(data = sim, archetype_formula = sam_form, species_formula = sp_form)
}

make_rcp_data <- function(family = "bernoulli", nRCP = 2, S = 6, n = 80,
                           seed = 1, power = 1.6) {
  set.seed(seed)
  x1 <- stats::runif(n, -2, 2)
  Xmat <- cbind(1, x1)
  form.RCP <- stats::as.formula(paste0("cbind(",
    paste0("spp", 1:S, collapse = ","), ")~x1"))
  sim <- suppressMessages(regional_mix.simulate(nRCP = nRCP, S = S, p.x = 2, p.w = 0, n = n,
    alpha = rnorm(S, 1, 0.3),
    tau = matrix(rnorm((nRCP - 1) * S, 0, 0.5), nrow = nRCP - 1),
    beta = matrix(rnorm((nRCP - 1) * 2, 0, 0.3), nrow = nRCP - 1),
    X = Xmat, family = family,
    logDisps = if (family %in% c("negative.binomial", "tweedie", "gaussian")) rep(0, S) else NULL,
    powers = if (family == "tweedie") rep(power, S) else NULL))
  sim$x1 <- x1
  list(data = sim, rcp_formula = form.RCP)
}

quiet_control <- list(quiet = TRUE)
