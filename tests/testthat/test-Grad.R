test_that("input validation", {
  expect_error(Grad(x = 1:4, FUN = "rubbish"), "must be a function")
  expect_error(Grad(x = 1:4, FUN = sum, h = c(0.01, 0.02)), "must have length")
  expect_error(Grad(x = 1:4, FUN = sum, side = c(0, 1, 2, -2)), "must contain values")
  expect_error(Grad(x = 1:4, FUN = sum, h = 0), "must be positive")
  expect_error(Grad(x = 1:4, FUN = sum, h = -0.001), "must be positive")
  expect_error(Grad(x = 1:4, FUN = function(x) "0.1"), "numeric values only")
  expect_error(Grad(x = 1:4, FUN = sum, side = c(-1, 1)), "'side' argument must")
  expect_error(Grad(x = 1:4), "Pass the function")
  expect_error(Grad(x = 1:4, FUN = sum, deriv.order = c(1, 2)), "'deriv.order' must have length")
  expect_error(Grad(x = 1:4, FUN = sum, acc.order = c(1, 2)), "'acc.order' must have length")
})

test_that("compatibility with numDeriv", {
  expect_warning(Grad(x = 1:4, func = sum), "Use the argument")
  expect_error(Grad(x = 1:4, func = sum, method = "complex"), "Complex derivatives not implemented")
  expect_equal(Grad(x = 1:4, FUN = sum, side = NULL, report = 0), rep(1, 4), tolerance = 1e-10)
})

test_that("Gradients are correct", {
  f <- function(x) sum(sin(x))
  expect_equal(Grad(x = 1:4, FUN = f, report = 0), cos(1:4), tolerance = 1e-10)
})

test_that("Jacobians are correct", {
  f <- function(x) c(sine = sum(sin(x)), expon = sum(exp(x)))
  expect_equal(Grad(x = 1:3, f, report = 0), rbind(sine = cos(1:3), expon = exp(1:3)),
               tolerance = 1e-9)

  x <- structure(1:3, names = c("A", "B", "C"))
  g <- Grad(f, x)
  expect_equal(colnames(g), names(x))
  expect_equal(rownames(g), c("sine", "expon"))
  expect_equal(attr(g, "step.size.method"), "default")
  expect_equal(attr(Grad(f, x, h = 0.01), "step.size.method"), "user-supplied")
})

test_that("Compatibility with numDeriv syntax", {
  f <- function(x) sum(sin(x))
  expect_equal(Grad(x = 1:3, FUN = f, method = "simple"), Grad(x = 1:3, f), tolerance = 1e-9)
})

test_that("Gradient step is auto-selected well", {
  f <- function(x) sum(sin(x))
  expect_equal(attr(Grad(x = 1:3, FUN = f), "step.size.method"), "default")
  expect_equal(attr(Grad(x = 1:3, FUN = f, h = 0.01), "step.size.method"), "user-supplied")
  expect_equal(attr(Grad(x = 1:3, FUN = f, h = "SW"), "step.size.method"), "SW")
  expect_equal(attr(Grad(x = 1:3, FUN = f, h = "CR"), "step.size.method"), "CR")

  g <- Grad(x = 1:3, FUN = f, h = "SW", report = 2)
  expect_equal(attr(g, "step.search")[["exitcode"]], rep(0, 3), tolerance = 1e-15)
  expect_length(attr(g, "step.search")[["iterations"]], 3)
})

test_that("parallelisation of Grad works", {
  expect_equal(Grad(x = 1:3, FUN = sum, cores = 1), Grad(x = 1:3, FUN = sum, cores = 2))
})

