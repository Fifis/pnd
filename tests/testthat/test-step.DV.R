test_that("Dumontet-Vignes step selection handles inputs well", {
  f <- function(x) return(NA)
  expect_error(step.DV(x = 2, f), "Could not compute the function value")
  expect_error(step.DV(sin, 1, range = c(0, 1)), "must be a positive vector of length 2")
})

test_that("Dumontet-Vignes step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.DV(x = 2, f)
  expect_identical(s$exitcode, 0L)
  expect_lt(sum(s$abs.error), 1e-6)
  expect_equal(s$value, 32, tolerance = 1e-8)
  u <- s$iterations$ratio[length(s$iterations$ratio)]
  u <- max(u, 1/u)
  # Stopping criterion
  expect_gte(u, 2)
  expect_lte(u, 15)

  s2 <- step.DV(x = 2, f, range = c(1e-10, 1e-7))
  expect_identical(s2$exitcode, 3L)

  s3 <- step.DV(x = 2, f, range = c(1e-3, 1e-1), maxit = 10)
  expect_identical(s3$exitcode, 6L)

  s4 <- step.DV(x = 2, f, h0 = 1000, maxit = 5)
  expect_identical(s4$exitcode, 6L)

  # Too large a size must be limited -- the range must be over-ridden
  s5 <- step.DV(x = 2, f, h0 = 1000, range = c(1e2, 1e4), maxit = 5)
  expect_equal(s5$par, 0.2, tolerance = 1e-12)
  expect_identical(s5$exitcode, 6L)

  expect_identical(step.DV(x = 2, f, maxit = 1)$exitcode, 7L)
})

test_that("Dumontet--Vignes algorithm stops if the function returns NA for all allowed step sizes", {
  f <- function(x) ifelse(abs(x - 2) < 1e-8, x^4, NA)
  expect_error(step.DV(f, 2, range = c(1e-7, 1e-2)), "attempts of step shrinkage")
})

test_that("Tweaking the DV algorithm for noisier functions", {
  f <- function(x) x^4
  s.perfect <- step.DV(x = 2, f, h0 = 1e-7, max.rel.error = 2e-16)
  s.noisy <- step.DV(x = 2, f, h0 = 1e-7, max.rel.error = 2e-8)
  expect_lt(s.perfect$par, s.noisy$par)
})

test_that("DV for functions with near-zero f''' stops immediately", {
  # Quadratic function, f''' = 0
  s1 <- step.DV(function(x) x^2, 1)
  expect_lte(s1$counts, 2)
  expect_identical(s1$exitcode, 1L)

  s2 <- step.DV(function(x) pi*x + exp(1), 1)
  expect_lte(s2$counts, 2)
  expect_identical(s2$exitcode, 1L)
})

test_that("Parallelisation in DV works", {
  expect_identical(step.DV(sin, 1, cores = 1), step.DV(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.DV(sin, 1, cores = 1), step.DV(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
