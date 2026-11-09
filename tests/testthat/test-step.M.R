test_that("Mathur's AutoDX handles inputs well", {
  expect_error(step.M(sin, 1, range = c(0, 1), cores = 1), "must be a positive vector of length 2")
  expect_warning(step.M(sin, 1, range = c(1e-4, 1e-6), cores = 1), "was extended to")
})

test_that("Mathur's step selection behaves reasonably", {
  s <- step.M(x = pi/4, sin, cores = 1)
  expect_identical(s$exitcode, 0L)
  expect_equal(s$value, sqrt(2)/2, tolerance = 1e-8)
})

test_that("Mathur's general-order estimates behave reasonably", {
  s <- suppressWarnings(step.M(function(z) z^6, 0, deriv.order = 2, acc.order = 4))
  h <- s$iterations$h
  error <- s$iterations$est.error[, "trunc"]
  expect_equal(mean(error/h^4, na.rm = TRUE), 8, tolerance = 1e-12)
  expect_true(all(s$iterations$aprime == 2L))
})

test_that("Mathur interpolates the estimated d^2f/dx^2 without new evaluations", {
  calls <- 0L
  f <- function(z) { calls <<- calls + 1L; z^4 }
  s <- step.M(f, 1, deriv.order = 2, cores = 1)
  expect_equal(s$value, 12 + 2*s$par^2, tolerance = 1e-12)
  expect_equal(calls, 2*s$counts + 1)
})

test_that("Mathur's fallback uses the correct rounding target", {
  f <- function(z) ifelse(z == 0 | abs(z) %in% c(1/4, 1/2, 1, 2), 1 + z^6, NA_real_)
  p <- 1e-8
  s <- suppressWarnings(step.M(f, 0, deriv.order = 2, acc.order = 4, max.rel.error = p))
  # K*|f^(6)|=8, |weights| sum to 16/3, and the local function maximum is 65/64
  B <- (16*p/3 + .Machine$double.eps/2) * 65/64
  expect_equal(s$iterations$fallback$target, (16*B^2)^(1/3), tolerance = 1e-12)
})

test_that("Mathur's step selection respects non-integer stencil offsets", {
  expect_equal(step.M(exp, 1, shrink.factor = 0.6, cores = 1)$value, exp(1), tolerance = 1e-7)
  expect_equal(step.M(exp, 1, shrink.factor = 0.8, cores = 1)$value, exp(1), tolerance = 1e-7)
  expect_equal(step.M(exp, 1, shrink.factor = 1/sqrt(2), cores = 1)$value, exp(1), tolerance = 1e-7)

  # Higher accuracy order also works
  expect_equal(suppressWarnings(step.M(identity, 1, acc.order = 4, shrink.factor = 0.6, cores = 1))$value, 1, tolerance = 1e-10)
})

test_that("Mathur's algorithm returns reasonable values even with bad slopes", {
  # TODO: come up with an example where the slopes is slightly off

  f <- function(x) ifelse(x %in% 1:2, x^2, NA)
  expect_warning(m <- step.M(f, 1, cores = 1), "<3 finite function values")
  expect_identical(m$exitcode, 3L)
})

test_that("Mathur's algorithm may return unreasonable values", {
  f <- function(x) x^3 + 1/x  # Produces wildly wrong results due to the shape of the error plot
  expect_gt(plot(step.M(f, 1, cores = 1))$par, 0.01)
})

test_that("Mathur's algorithm returns expected non-zero exit codes", {
  # Noisy right branch
  expect_identical(suppressWarnings(step.M(function(x) x^2, x = 1e-8)$exitcode), 2L)
  # No branches at all
  expect_identical(suppressWarnings(step.M(function(x) x^2, x = 0)$exitcode), 2L)
})

test_that("Parallelisation in Mathur's algorithm works", {
  expect_identical(step.M(sin, 1, cores = 1), step.M(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.M(sin, 1, cores = 1), step.M(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
