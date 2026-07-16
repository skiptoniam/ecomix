test_that("species_mix fits bernoulli data and returns a sane object", {
  d <- make_sam_data("bernoulli")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "bernoulli", nArchetypes = 2, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.finite(fm$logl))
  expect_true(is.list(fm$coefs))
  expect_length(fm$alpha, 6)
  expect_equal(dim(fm$coefs$beta), c(2, 2))
  expect_equal(sum(fm$pi), 1, tolerance = 1e-6)
})

test_that("species_mix fits poisson data", {
  d <- make_sam_data("poisson")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "poisson", nArchetypes = 2, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.finite(fm$logl))
})

test_that("species_mix fits negative.binomial data and optimises dispersion", {
  d <- make_sam_data("negative.binomial")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "negative.binomial", nArchetypes = 2, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.finite(fm$logl))
  expect_true(all(is.finite(fm$theta)))
  expect_true(all(fm$theta > 0))
})

test_that("species_mix fits gaussian data", {
  d <- make_sam_data("gaussian")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "gaussian", nArchetypes = 2, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.finite(fm$logl))
  expect_true(all(fm$theta > 0))
})

test_that("species_mix fits tweedie data", {
  d <- make_sam_data("tweedie")
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "tweedie", power = 1.6, nArchetypes = 2, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.finite(fm$logl))
  expect_true(all(fm$theta > 0))
})

test_that("species_mix works for a single archetype (G=1)", {
  d <- make_sam_data("negative.binomial", nArchetypes = 1)
  fm <- species_mix(archetype_formula = d$archetype_formula, species_formula = d$species_formula,
                     data = d$data, family = "negative.binomial", nArchetypes = 1, control = quiet_control)
  expect_s3_class(fm, "species_mix")
  expect_true(is.list(fm$coefs))
  expect_equal(unname(fm$pi), 1)
  expect_true(all(fm$theta > 0))
  # theta must actually be optimised, not frozen at the EM starting value
  expect_false(isTRUE(all.equal(as.numeric(fm$theta), rep(fm$theta[1], length(fm$theta)))))
})

test_that("species_mix.fit errors informatively without an archetype_formula", {
  d <- make_sam_data("bernoulli")
  expect_message(fm <- species_mix(NULL, d$species_formula, data = d$data,
                                    family = "bernoulli", nArchetypes = 2))
  expect_null(fm)
})
