test_that("Stepleman-Winarsky step selection handles inputs well", {
  f <- function(x) return(NA)
  expect_error(step.SW(x = 2, f), "must be finite")

  f2 <- function(x) if (x == 2) return(3) else return(NA)
  expect_error(step.SW(x = 2, f2), "is finite -- reduce")
})


test_that("Stepleman-Winarsky step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.SW(x = 2, f, diagnostics = TRUE)
  expect_equal(s$exitcode, 0)
  expect_lt(s$abs.error, 1e-6)
  expect_equal(s$value, 32, tolerance = 1e-8)
  monot <- s$iterations$monotone[sum(s$counts), ]
  # Stopping criterion
  expect_true(any(!monot))
})

test_that("SW algorithm detects if h0 is too low", {
  f <- function(x) x^4
  s <- step.SW(x = 2, f, h0 = 1e-9, diagnostics = TRUE)
  expect_equal(s$exitcode, 0)
  expect_lt(s$iterations$h[1], s$iterations$h[2])
  expect_equal(s$value, 32, tolerance = 1e-8)
})

test_that("SW algorithm detects if h0 is too high", {
  f <- function(x) x^4
  s <- step.SW(x = 2, f, h0 = 10, diagnostics = TRUE)
  expect_equal(s$exitcode, 0)
  expect_gt(s$counts["preliminary"], 3)

  s2 <- step.SW(x = 2, f, h0 = 10, shrink.factor = 4, diagnostics = TRUE)
  expect_equal(s2$counts["preliminary"], s$counts["preliminary"])
  expect_lt(s2$counts["main"], s$counts["main"])
})

test_that("SW fails when a large h0 invalidates the est. trunc. error", {
  expect_warning(step.SW(x = pi/4, sin, h0 = 1000, diagnostics = TRUE),
                 "exceeds the absolute value")
})

