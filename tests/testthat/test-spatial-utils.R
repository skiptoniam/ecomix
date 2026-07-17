test_that("covariance_function returns sig2 at distance 0 and decays with distance for each nu", {
  for (nu in c(1/2, 3/2, 5/2)) {
    cov0 <- ecomix:::covariance_function(0, sig2 = 2, rho = 0.5, nu = nu)
    cov1 <- ecomix:::covariance_function(1, sig2 = 2, rho = 0.5, nu = nu)
    expect_equal(cov0, 2)
    expect_lt(cov1, cov0)
  }
})

test_that("covariance_function errors for an unsupported nu", {
  expect_error(ecomix:::covariance_function(1, sig2 = 1, rho = 0.5, nu = 1),
               "nu must be one one of")
})

test_that("extend pads a vector to a power-of-2 length at least `factor` times as long", {
  x <- seq(0, 10, by = 1)
  e <- ecomix:::extend(x, factor = 2)
  expect_equal(length(e), 32)
  expect_true(length(e) >= 2 * length(x))
  expect_equal(log2(length(e)) %% 1, 0)
  idx <- attr(e, "idx")
  expect_equal(e[idx[1]:idx[2]], x)
})

test_that("seq_range reproduces seq() over a two-element range", {
  expect_equal(ecomix:::seq_range(c(1, 5)), seq(1, 5, by = 1))
  expect_equal(ecomix:::seq_range(c(0, 1), by = 0.25), seq(0, 1, by = 0.25))
})

test_that("block_circulant_basis returns one distance-based value per grid cell", {
  x <- y <- seq(-1, 1, len = 4)
  bcb <- ecomix:::block_circulant_basis(x, y, sig2 = 1, rho = 0.5)
  expect_length(bcb, length(x) * length(y))
  expect_true(all(is.finite(bcb)))
  expect_equal(bcb[1], max(bcb))
})

test_that("rGMRF simulates a finite field of the requested dimensions, with and without a nugget", {
  set.seed(1)
  x <- y <- seq(-2, 2, len = 16)
  sim1 <- rGMRF(x, y, sig2 = 1, rho = 0.5, nu = 1/2)
  expect_equal(dim(sim1), c(length(x), length(y)))
  expect_true(all(is.finite(sim1)))

  sim2 <- rGMRF(x, y, sig2 = 1, rho = 0.5, nu = 1/2, nugget = 0.1)
  expect_equal(dim(sim2), dim(sim1))
  expect_true(all(is.finite(sim2)))
})

test_that("rGMRF preserves NAs through the FFT round trip", {
  set.seed(1)
  fs <- ecomix:::setupFFTgp(x = seq(-1, 1, len = 4), y = seq(-1, 1, len = 4))
  meanmat <- matrix(rnorm(16), 4, 4)
  meanmat[2, 2] <- NA
  out <- ecomix:::fft_GP(meanmat, fs)
  expect_true(is.na(out[2, 2]))
  expect_true(all(is.finite(out[-c(2), -c(2)])))
})
