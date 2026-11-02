test_that("Mathur's AutoDX handles inputs well", {
  expect_error(step.M(sin, 1, range = c(0, 1), cores = 1), "must be a positive vector of length 2")
  expect_warning(step.M(sin, 1, range = c(1e-4, 1e-6), cores = 1), "was extended to")
})

test_that("Mathur's step selection behaves reasonably", {
  s <- step.M(x = pi/4, sin, cores = 1)
  expect_identical(s$exitcode, 0L)
  expect_equal(s$value, sqrt(2)/2, tolerance = 1e-8)
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
