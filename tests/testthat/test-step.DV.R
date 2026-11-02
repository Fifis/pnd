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

test_that("DV exit messages match the termination status", {
  single <- step.DV(sin, 1, maxit = 1)
  expect_identical(single$exitcode, 7L)

  capped <- step.DV(function(x) x^3, 1, h0 = 1, range = c(1e-10, 1), max.rel.error = 0.1, maxit = 1)
  expect_identical(capped$exitcode, 4L)

  boundary <- step.DV(function(x) x^4, 2, range = c(1e-3, 1e-1), maxit = 10)
  expect_identical(boundary$exitcode, 6L)
  expect_match(boundary$message, "left end")

  right <- step.DV(function(x) x^3, 1, h0 = 1e-3, range = c(1e-5, 1e-3), max.rel.error = 0.01, maxit = 2)
  expect_identical(right$exitcode, 6L)
  expect_match(right$message, "right end")
})

test_that("DV uses the last evaluated spacing when iterations are exhausted", {
  scale <- .Machine$double.eps^(-2/15)
  s <- step.DV(function(x) x^3, 1, h0 = (1/8)/scale, range = c(1/128, 8)/scale, max.rel.error = 1e-8, maxit = 2)
  k <- tail(s$iterations$k, 1)
  # For x^3 at x = 1, f''' = 6 and the mean of f(1-k), f(1+k) is 1 + 3*k^2.
  expected <- (1.68 * 1e-8 * (1 + 3*k^2) / 6)^(1/3)
  expect_identical(s$counts, 2L)
  expect_identical(s$exitcode, 5L)
  expect_equal(s$par, expected, tolerance = 1e-8)
})

test_that("DV retains convergence on the last permitted iteration", {
  limits <- c(1.5, 10)
  s <- step.DV(function(x) x^3, 1, h0 = 1e-3, range = c(1e-5, 0.1), max.rel.error = 1e-6, ratio.limits = limits, maxit = 2)
  expect_identical(s$counts, 2L)
  expect_gte(tail(s$iterations$ratio, 1), limits[1])
  expect_lte(tail(s$iterations$ratio, 1), limits[2])
  expect_identical(s$exitcode, 0L)
})

test_that("DV history records accepted evaluations after shrinking", {
  f <- function(x) ifelse(abs(x-1) < 0.5, x^3, NA_real_)
  s <- step.DV(f, 1, h0 = 0.01, maxit = 1)
  expect_true(all(is.finite(s$iterations$f)))
  expect_true(all(abs(s$iterations$x-1) < 0.5))
  expect_equal(unname(s$iterations$f[1, ]), unname(s$iterations$x[1, ]^3), tolerance = 1e-10)
  expect_equal(unname(s$iterations$deriv3[1, "f3"]), 6, tolerance = 1e-8)
})
