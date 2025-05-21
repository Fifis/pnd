test_that("Kostyrka's step selection handles inputs well", {
  expect_error(step.K(sin, 1, range = c(0, 1), cores = 1), "must be a positive vector of length 2")
  expect_warning(step.K(sin, 1, range = c(1e-4, 1e-6), cores = 1), "was extended to")
})

test_that("Kostyrka's step selection behaves reasonably", {
  s <- step.K(x = pi/4, sin, cores = 1)
  expect_identical(s$exitcode, 0L)
  expect_equal(s$value, sqrt(2)/2, tolerance = 1e-10)
})

test_that("Kostyrka's step selection works for gradients", {
  f <- function(x) sin(x[1]) * cos(x[2]) + exp(x[3])
  s <- gradstep(f, 1:3, method = "K")
  expect_identical(s$exitcode, rep(0L, 3))
})

test_that("Kostyrka's method returns reasonable values with unfortunate inputs", {
  s <- step.K(sin, 1, shrink.factor = 1/8, max.rel.error = 0.5, cores = 1)
  expect_equal(s$value, cos(1), tolerance = 1e-9)
  expect_identical(s$exitcode, 0L)

  f <- function(x) ifelse(x %in% 0:2, x^2, NA)
  expect_identical(step.K(f, 1, cores = 1)$exitcode, 2L)
})

test_that("Kostyrka's method returns reasonable values with tricky functions", {
  s <- step.K(function(x) x^3+1/x, 1)
  expect_identical(s$exitcode, 0L)
  s <- step.K(function(x) 6*x^5 - 56*x^3 + 1/x^2, 1)
  expect_identical(s$exitcode, 0L)
  expect_equal(s$value, -140, tolerance = 1e-10)
})

test_that("Mathur's algorithm returns expected non-zero exit codes", {
  # Very noisy right branch
  expect_identical(suppressWarnings(step.K(function(x) x^2, x = 1e-8)$exitcode), 1L)
  # No right branch at all
  expect_identical(suppressWarnings(step.K(function(x) x^2, x = 0)$exitcode), 2L)
})

test_that("Parallelisation in Kostyrka's algorithm works", {
  expect_identical(step.K(sin, 1, cores = 1), step.K(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.K(sin, 1, cores = 1), step.K(sin, 1, cl = clus))
  parallel::stopCluster(clus)

  # Testing a slow function
  # f <- function(x) {Sys.sleep(0.1); sin(x)}
  # system.time(step.K(f, 1))
  # system.time(step.K(f, 1, cores = 12))
})
