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
