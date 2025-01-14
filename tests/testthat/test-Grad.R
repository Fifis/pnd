test_that("compatibility with numDeriv", {
  expect_warning(Grad(x = 1:4, func = sum, vectorised = FALSE), "Use the argument")
  expect_warning(Grad(x = 1:3, FUN = sum, method = "simple", vectorised = FALSE), "numDeriv-like syntax")
  expect_equal(suppressWarnings(Grad(x = 1:3, FUN = sum, method = "simple", vectorised = FALSE, report = 0)),
               Grad(x = 1:3, sum, h = 1:3 * 2 * sqrt(.Machine$double.eps), vectorised = FALSE, report = 0),
               tolerance = 1e-15)
  expect_error(Grad(x = 1:4, func = sum, method = "complex", vectorised = FALSE), "Complex derivatives not implemented")
  expect_equal(Grad(x = 1:4, FUN = sum, side = NULL, vectorised = FALSE, report = 0), rep(1, 4), tolerance = 1e-10)
  # TODO: Richardson
})

test_that("gradients are correct", {
  f <- function(x) sum(sin(x))
  expect_equal(Grad(x = 1:4, f, vectorised = FALSE, report = 0), cos(1:4), tolerance = 1e-10)

  x <- structure(1:3, names = c("A", "B", "C"))
  g <- Grad(f, x, h = 0.01, vectorised = FALSE)
  expect_equal(attr(g, "step.size.method"), "user-supplied")
  expect_equal(names(g), c("A", "B", "C"))
})

test_that("Gradient step is auto-selected well", {
  f <- function(x) sum(sin(x))
  expect_equal(attr(Grad(x = 1:3, f, vectorised = FALSE), "step.size.method"), "default")
  expect_equal(attr(Grad(x = 1:3, f, h = 0.01, vectorised = FALSE), "step.size.method"), "user-supplied")
  expect_equal(attr(Grad(x = 1:3, f, h = "SW", vectorised = FALSE), "step.size.method"), "SW")
  expect_equal(attr(Grad(x = 1:3, f, h = "CR", vectorised = FALSE), "step.size.method"), "CR")

  g <- Grad(x = 1:3, FUN = f, h = "SW", vectorised = FALSE, report = 2)
  expect_equal(attr(g, "step.search")[["exitcode"]], rep(0, 3), tolerance = 1e-15)
  expect_length(attr(g, "step.search")[["iterations"]], 3)
})

test_that("parallelisation of Grad works", {
  expect_equal(Grad(x = 1:3, FUN = sum, vectorised = FALSE, cores = 1),
               Grad(x = 1:3, FUN = sum, vectorised = FALSE, cores = 2))
})

test_that("function dimension check works", {
  f <- function(x) c(sum(sin(x)), sum(exp(x)))
  expect_error(Grad(f, 1:3, vectorised = FALSE), "vector-valued function")
  expect_error(Grad(f, 1:3, h = "SW", vectorised = FALSE), "returns a scalar")

  w <- expect_error(pnd::Grad(sin, x = 1:4, deriv.order = 1, vectorised = FALSE))
})
