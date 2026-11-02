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

test_that("SW retains zero differences during an otherwise defined search", {
  s <- step.SW(exp, 1, h0 = 1e-20, range = c(1e-24, 0.1), cores = 1)
  expect_identical(s$iterations$deriv[1], 0)
  expect_equal(s$value, exp(1), tolerance = 1e-8)
  expect_identical(s$exitcode, 0L)
  expect_identical(step.SW(function(z) 1+z^3, 0, cores = 1)$value, 0)
})

test_that("SW reports an undefined cancellation measure at zero", {
  msg <- "cancellation measure is undefined"
  expect_error(step.SW(function(z) z^2, 0, cores = 1), msg)
  expect_error(step.SW(function(z) (z-2)^2, 2, h0 = 1/4, cores = 1), msg)
})

test_that("SW follows the paper's shrinking branch for zero endpoints", {
  expect_identical(step.SW(function(z) 0, 1, cores = 1)$value, 0)
  s <- step.SW(function(z) z+z^2, 0, h0 = 1, cores = 1)
  expect_equal(s$value, 1)
  expect_identical(s$exitcode, 0L)
  # Equal sampled values do not imply a zero derivative: here f'(0) = -1.
  s <- step.SW(function(z) z*(z^2-1), 0, h0 = 1, cores = 1)
  expect_identical(s$iterations$deriv[1], 0)
  expect_equal(s$value, -1, tolerance = 1e-8)
})

test_that("SW initial recovery evaluates and records the reduced spacing", {
  seen <- numeric()
  f <- function(z) {
    seen <<- c(seen, z)
    if (abs(z) <= 1/8) exp(z) else NA_real_
  }
  s <- step.SW(f, 0, h0 = 1/2, range = c(1e-12, 1), cores = 1)
  expect_equal(seen[1:9], c(-1/2, 0, 1/2, -1/4, 0, 1/4, -1/8, 0, 1/8))
  expect_equal(s$iterations$h[1:3], c(1/8, 1/16, 1/32))
  expect_equal(s$iterations$x[1, ], c(-1/8, 1/8))
  expect_equal(s$iterations$f[1, ], exp(c(-1/8, 1/8)))
  expect_equal(s$iterations$deriv[1], sinh(1/8)/(1/8), tolerance = 1e-10)
  etrunc <- abs((sinh(1/16)/(1/16) - sinh(1/8)/(1/8)) / (1 - 2^2))
  expect_equal(unname(s$iterations$est.error[2, "trunc"]), etrunc, tolerance = 1e-10)
  expect_equal(s$value, 1, tolerance = 1e-8)
})

test_that("SW initial recovery tries the lower boundary before stopping", {
  seen <- numeric()
  f <- function(z) {
    seen <<- c(seen, z)
    if (z == 1) 1 else NA_real_
  }
  expect_error(step.SW(f, 1, h0 = 1/2, range = c(1/8, 1), cores = 1), "attempts of step shrinkage")
  expect_equal(seen, c(1/2, 1, 3/2, 3/4, 1, 5/4, 7/8, 1, 9/8))
})


test_that("SW fails when a large h0 invalidates the est. trunc. error", {
  expect_identical(step.SW(x = pi/4, FUN = sin, h0 = 1000)$exitcode, 4L)
})

test_that("SW capped results use the evaluated final step", {
  x <- pi/4
  for (limit in c(40L, 4L)) {
    s <- step.SW(sin, x, h0 = 1000, maxit = limit, cores = 1)
    expect_identical(s$exitcode, if (limit == 4L) 5L else 4L)
    expect_equal(s$value, cos(x)*sin(s$par)/s$par, tolerance = 1e-10)
    last <- length(s$iterations$h)
    expect_identical(s$par, s$iterations$h[last])
    expect_identical(s$value, s$iterations$deriv[last])
    expect_equal(unname(s$abs.error), unname(s$iterations$est.error[last, ]), tolerance = 1e-10)
  }
})

test_that("SW retains its upper-bound result when the iteration limit is reached", {
  d10 <- sin(10)/10
  d5 <- sin(5)/5
  expected.error <- c(trunc = abs((d5-d10) / (1-2^2)) / 0.5,  round = 0.5*abs(sin(10))/10)
  for (limit in c(40L, 3L)) {
    s <- suppressWarnings(step.SW(sin, 0, h0 = 10, range = c(0.001, 10), max.rel.error = 0.5, maxit = limit, cores = 1))
    expect_identical(s$exitcode, if (limit == 3L) 5L else 3L)
    expect_equal(s$par, 10)
    expect_equal(s$value, d10, tolerance = 1e-10)
    expect_equal(unname(s$abs.error), unname(expected.error), tolerance = 1e-10)
  }
})

test_that("Parallelisation in SW works", {
  expect_identical(step.SW(sin, 1, cores = 1), step.SW(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.SW(sin, 1, cores = 1), step.SW(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
