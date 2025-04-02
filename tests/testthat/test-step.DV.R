test_that("Dumontet-Vignes step selection handles inputs well", {
  f <- function(x) return(NA)
  expect_error(step.DV(x = 2, f), "Could not compute the function value")
  expect_error(step.DV(sin, 1, range = c(0, 1)), "must be a positive vector of length 2")
})

test_that("Dumontet-Vignes step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.DV(x = 2, f, diagnostics = TRUE)
  expect_equal(s$exitcode, 0)
  expect_lt(s$abs.error, 1e-6)
  expect_equal(s$value, 32, tolerance = 1e-8)
  u <- s$iterations$ratio[length(s$iterations$ratio)]
  u <- max(u, 1/u)
  # Stopping criterion
  expect_gte(u, 2)
  expect_lte(u, 15)

  s2 <- step.DV(x = 2, f, range = c(1e-10, 1e-7))
  expect_equal(s2$exitcode, 3)
  expect_true(grepl("too close to the right end", s2$message))

  s3 <- step.DV(x = 2, f, range = c(1e-3, 1e-1))
  expect_equal(s3$exitcode, 5)
  expect_true(grepl("on the right end", s3$message))

  s4 <- step.DV(x = 2, f, h0 = 100, maxit = 5)
  expect_equal(s4$exitcode, 5)

  expect_equal(step.DV(x = 2, f, maxit = 1)$exitcode, 6)
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
  expect_equal(s1$exitcode, 1)

  s2 <- step.DV(function(x) pi*x + exp(1), 1)
  expect_lte(s2$counts, 2)
  expect_equal(s2$exitcode, 1)
})

test_that("Parallelisation in DV works", {
  expect_equal(step.DV(sin, 1, cores = 1), step.DV(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_equal(step.DV(sin, 1, cores = 1), step.DV(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
