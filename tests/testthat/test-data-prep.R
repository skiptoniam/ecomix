test_that("table_to_species_data aggregates a long-format table into a sites x species matrix", {
  obs <- data.frame(
    site_id = c("s1", "s1", "s2", "s2", "s3"),
    species_id = c("a", "b", "a", "c", "b"),
    count = c(3, 1, 2, 5, 4)
  )
  m <- table_to_species_data(obs, site_id = "site_id", species_id = "species_id",
                              measurement_id = "count")
  expect_equal(dim(m), c(3, 3))
  expect_equal(rownames(m), c("s1", "s2", "s3"))
  expect_equal(colnames(m), c("a", "b", "c"))
  expect_equal(m["s1", "a"], 3)
  expect_equal(m["s2", "c"], 5)
  expect_equal(m["s1", "c"], 0)
})

test_that("table_to_species_data treats missing measurement_id as presence/absence", {
  obs <- data.frame(
    site_id = c("s1", "s1", "s2"),
    species_id = c("a", "b", "a")
  )
  m <- table_to_species_data(obs, site_id = "site_id", species_id = "species_id")
  expect_equal(m["s1", "a"], 1)
  expect_equal(m["s1", "b"], 1)
  expect_equal(m["s2", "b"], 0)
})

test_that("table_to_species_data drops sites with blank or NA ids and all-zero rows", {
  obs <- data.frame(
    site_id = c("s1", NA, "s2", ""),
    species_id = c("a", "a", "a", "a"),
    count = c(1, 99, 2, 99)
  )
  m <- table_to_species_data(obs, site_id = "site_id", species_id = "species_id",
                              measurement_id = "count")
  expect_setequal(rownames(m), c("s1", "s2"))
})

test_that("make_mixture_data combines species and covariate data with an intercept column", {
  sp <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, dimnames = list(NULL, c("spA", "spB")))
  cov <- data.frame(x1 = c(0.1, 0.2, 0.3))
  out <- make_mixture_data(sp, cov)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("spA", "spB", "const", "x1"))
  expect_equal(out$const, rep(1, 3))
  expect_equal(nrow(out), 3)
})

test_that("make_mixture_data errors when species and covariate row counts differ", {
  sp <- matrix(c(1, 0, 1, 0), nrow = 2, dimnames = list(NULL, c("spA", "spB")))
  cov <- data.frame(x1 = c(0.1, 0.2, 0.3))
  expect_error(make_mixture_data(sp, cov), "dimensions of species matrix")
})

test_that("make_mixture_data rejects non-finite or wrong-typed inputs", {
  sp_bad <- matrix(c(1, NA, 1, 0), nrow = 2, dimnames = list(NULL, c("spA", "spB")))
  cov <- data.frame(x1 = c(0.1, 0.2))
  expect_error(make_mixture_data(sp_bad, cov))
  expect_error(make_mixture_data(list(a = 1), cov))
})
