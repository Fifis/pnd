test_that("Stepleman--Winarsky step selection handles inputs well", {
  expect_error(step.SW(sin, 1, range = c(0, 1)), "must be a positive vector of length 2")
  f <- function(x) return(NA)
  expect_error(step.SW(x = 2, f), "must be finite")

  f2 <- function(x) if (x == 2) return(3) else return(NA)
  expect_error(step.SW(x = 2, f2), "attempts of step shrinkage")
})

test_that("Stepleman--Winarsky step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.SW(x = 2, f)
  expect_identical(s$exitcode, 0L)
  expect_lt(sum(s$abs.error), 1e-6)
  expect_equal(s$value, 32, tolerance = 1e-8)
  monot <- s$iterations$monotone[sum(s$counts), ]
  # Stopping criterion
  expect_false(all(monot))
})

test_that("SW algorithm detects if h0 is too low", {
  f <- function(x) x^4
  s <- step.SW(x = 2, f, h0 = 1e-9)
  expect_identical(s$exitcode, 0L)
  expect_lt(s$iterations$h[1], s$iterations$h[2])
  expect_equal(s$value, 32, tolerance = 1e-8)
})

test_that("SW algorithm takes longer if h0 is too high", {
  f <- function(x) x^4
  s <- step.SW(x = pi/2, f, h0 = 10)
  expect_identical(s$exitcode, 0L)
  expect_gt(s$counts["preliminary"], 3)

  s2 <- step.SW(x = pi/2, f, h0 = 1e4)
  expect_lt(s$counts["preliminary"], s2$counts["preliminary"])
})

test_that("SW algorithm for bad ranges", {
  f <- function(x) x^4
  expect_identical(step.SW(x = 2, f, h0 = 0.1, range = c(0.01, 1))$exitcode, 3L)
  expect_identical(suppressWarnings(step.SW(x = 2, f, h0 = 1e-9, range = c(1e-10, 1e-8)))$exitcode, 3L)
})

test_that("SW algorithm stops if the function returns NA for all allowed step sizes", {
  f <- function(x) ifelse(abs(x - 2) < 1e-8, x^4, NA)
  expect_error(step.SW(f, 2, range = c(1e-7, 1e-2)), "attempts of step shrinkage")
})


test_that("SW fails when a large h0 invalidates the est. trunc. error", {
  expect_identical(step.SW(x = pi/4, FUN = sin, h0 = 1000)$exitcode, 4L)
})

test_that("Parallelisation in SW works", {
  expect_identical(step.SW(sin, 1, cores = 1), step.SW(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.SW(sin, 1, cores = 1), step.SW(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
