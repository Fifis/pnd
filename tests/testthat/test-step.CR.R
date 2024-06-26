test_that("Curtis-Reid step selection handles inputs well", {
  f <- function(x) return(NA)
  expect_error(step.CR(x = 2, f), "Could not compute the function value")
  expect_error(step.CR(sin, 1, version = "orig"), "must be either")
  expect_error(step.CR(sin, 1, tol = 1e-4), "must be a positive number greater than 1")
  expect_warning(step.CR(sin, 1, acc.order = 4), "Setting acc.order")
  expect_equal(step.CR(sin, 1, range = c(0, 1))$exitcode, 0)
})


test_that("Curtis-Reid step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.CR(x = 2, f, diagnostics = TRUE)
  expect_equal(s$exitcode, 0)
  expect_lt(s$abs.error, 1e-5)
  expect_equal(s$value, 32, tolerance = 1e-8)
  u <- s$iterations$ratio[length(s$iterations$ratio)]
  expect_gt(u, 10)
  expect_lt(u, 1000)

  s2 <-  step.CR(x = 2, f, version = "modified")
  expect_lt(s2$abs.error, 1e-7)
  expect_equal(s2$value, 32, tolerance = 1e-8)

  s3 <- step.CR(x = 2, f, version = "modified", acc.order = 4)
  expect_lt(s2$abs.error, 5e-8)
  expect_equal(s3$value, 32, tolerance = 1e-8)

  s4 <- step.CR(x = sqrt(2), FUN = function(x) x^6 -2*x^4 - 4*x^2, h0 = 2^-16)
  expect_lt(abs(s4$value), 1e-8)

})

test_that("Curtis-Reid steps grow for linear functions", {
  expect_equal(step.CR(x = 1, function(x) x)$exitcode, 1)
})

test_that("Large and small initial values in Curtis--Reid cause range problems", {
  f <- function(x) x^4
  expect_equal(step.CR(x = 2, f, h0 = 1000)$exitcode, 3)
  expect_equal(step.CR(x = 2, f, h0 = 1000, range = c(1e-10, 1e5))$exitcode, 0)
  expect_equal(step.CR(x = 2, f, h0 = 1e-12)$exitcode, 1)  # In theory, it should be 3,
  # but the truncation error estimate is really zero...
  expect_equal(step.CR(x = 2, f, h0 = 9e-9, range = c(1e-8, 2e-8))$exitcode, 3)
  expect_equal(step.CR(x = 2, FUN = f, h0 = 1000, maxit = 4, range = c(1e-12, 1))$exitcode, 4)
})
