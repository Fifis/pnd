test_that("step size length h0 must be 1 or length(x)", {
  expect_error(gradstep(x = 1, FUN = sin, h0 = c(1e-5, 1e-6)), "must be a scalar")
  expect_error(gradstep(x = 1:2,
                        FUN = function(z) {if (length(z)>1) stop("Non-vectorised"); z^2}),
               "must be finite")
  expect_equal(gradstep(sin, 1, h0 = 1)$exitcode, 0)
  expect_equal(gradstep(function(x) sum(sin(x)), 1:3, h0 = 0.01)$exitcode, rep(0, 3))
})

test_that("step selection supports only scalar-valued functions", {
  expect_error(gradstep(x = 1, FUN = function(x) c(x, x^2)), "returns a scalar")
})

test_that("unsupported arguments withcause an error", {
  expect_error(gradstep(x = 1, FUN = sin, method = "CR", control = list(rubbish = TRUE)),
               "arguments are not supported")
  expect_error(gradstep(x = 1, FUN = sin, method = "CR", method.args = list(rubbish = TRUE)),
               "step-selection method arguments")
})

test_that("the range is correctly reversed", {
  expect_equal(gradstep(sin, 1, method = "CR", control = list(range = c(1e-6, 1e-8))),
               gradstep(sin, 1, method = "CR", control = list(range = c(1e-8, 1e-6))))
  expect_equal(gradstep(sin, 1, method = "DV", control = list(range = c(1e-6, 1e-8))),
               gradstep(sin, 1, method = "DV", control = list(range = c(1e-8, 1e-6))))
  expect_equal(gradstep(sin, 1, method = "SW", control = list(range = c(1e-6, 1e-8))),
               gradstep(sin, 1, method = "SW", control = list(range = c(1e-8, 1e-6))))
})

test_that("for unfortunate inputs, the search may hit the boundary", {
  s.sw <- suppressWarnings(gradstep(x = 1e10, FUN = sin, h0 = 1e-20, method = "SW"))
  expect_equal(s.sw$exitcode, 3)
  expect_true(is.finite(s.sw$abs.error))

  s.dv <- gradstep(x = 1, FUN = sin, method = "DV", control = list(range = c(1e-20, 1e-22)))
  expect_equal(s.dv$exitcode, 3)
})

test_that("gradstep correctly handles vector inputs", {
  f <- function(x) sum(sin(x))
  s.grad <- gradstep(x = 1:4, FUN = f,
                     control = list(diagnostics = TRUE))
  expect_equal(s.grad$exitcode, rep(0, 4))

  f1 <- function(x) sin(x) + sum(sin(2:4))
  s1 <- step.SW(f1, 1, diagnostics = TRUE)
  r <- s.grad$par[1] / s1$par
  r <- max(r, 1/r)
  expect_lte(r, 4)

  expect_error(gradstep(x = 1:4, f, h0 = c(1, 1)), "must have length 1 or length")
})

test_that("the error in vector inputs does not propagate too strongly", {
  s.grad <- gradstep(x = 1:4, FUN = function(x) sum(sin(x)),
                     control = list(diagnostics = TRUE))
  expect_equal(s.grad$exitcode, rep(0, 4))

  f1 <- function(x) sin(x) + sum(sin(2:4))
  s1 <- step.SW(f1, 1, diagnostics = TRUE)
  # Due to numerical errors, the number of steps can be different
  r <- s.grad$par[1] / s1$par
  r <- max(r, 1/r)
  expect_lte(r, 4)
})

test_that("gradstep accepts methods argument as a list", {
  expect_error(gradstep(x = 1:4, FUN = function(x) sum(sin(x)), diagnostics = TRUE))
  gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "SW",
           control = list())
})

test_that("gradstep correctly handles other inputs", {
  f <- function(x) sum(sin(x))
  expect_equal(gradstep(x = 1:4, FUN = f, h0 = NULL)$exitcode, rep(0, 4))
})

test_that("gradstep checks for conflicting arguments of FUN and methods", {
  # Root solver for sin(x) = a
  f <- function(x, tol) uniroot(function(y) sin(y) - x, c(0, 2), tol = tol)$root
  expect_error(gradstep(f, 0.5, tol = 1e-4, method = "CR"), "coincide with the arguments")
})
