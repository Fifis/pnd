test_that("Kostyrka's step selection handles inputs well", {
  expect_error(step.K(sin, 1, range = c(0, 1), cores = 1), "must be a positive vector of length 2")
  expect_warning(step.K(sin, 1, range = c(1e-4, 1e-6), cores = 1), "was extended to")
})

test_that("Kostyrka's step selection behaves reasonably", {
  s <- step.K(x = pi/4, sin, plot = TRUE, cores = 1)
  expect_equal(s$exitcode, 0)
  expect_equal(s$value, sqrt(2)/2, tolerance = 1e-10)
  if (file.exists("Rplot.pdf")) unlink("Rplot.pdf")
})

test_that("Kostyrka's method returns reasonable values with unfortunate inputs", {
  s <- step.K(sin, 1, shrink.factor = 1/8, max.rel.error = 0.5, cores = 1)
  expect_equal(s$value, cos(1), tolerance = 1e-9)
  expect_equal(s$exitcode, 0)

  f <- function(x) ifelse(x %in% 0:2, x^2, NA)
  expect_equal(step.K(f, 1, cores = 1)$exitcode, 2)
})

test_that("Kostyrka's method returns reasonable values with tricky functions", {
  s <- step.K(function(x) x^3+1/x, 1, plot = TRUE)
  expect_equal(s$exitcode, 0)
  s <- step.K(function(x) 6*x^5 - 56*x^3 + 1/x^2, 1, plot = TRUE)
  expect_equal(s$exitcode, 0)
  expect_equal(s$value, -140, tolerance = 1e-10)
})

test_that("Parallelisation in Kostyrka's algorithm works", {
  expect_equal(step.K(sin, 1, cores = 1), step.K(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_equal(step.K(sin, 1, cores = 1), step.K(sin, 1, cl = clus))
  parallel::stopCluster(clus)

  # Testing a slow function
  # f <- function(x) {Sys.sleep(0.1); sin(x)}
  # system.time(step.K(f, 1))
  # system.time(step.K(f, 1, cores = 12))
})
