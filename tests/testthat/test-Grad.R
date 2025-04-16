test_that("compatibility with numDeriv", {
  w <- testthat::capture_warnings(Grad(x = 1:4, func = sum))
  expect_match(w[1], "Use the argument")
  expect_match(w[2], "Use the argument")
  expect_warning(Grad(x = 1:3, FUN = sum, method = "simple"), "numDeriv-like syntax")
  s1 <- suppressWarnings(Grad(x = 1:3, FUN = sum, method = "simple"))
  s2 <- Grad(x = 1:3, sum, side = 1, acc.order = 1, h = 1e-5 * .Machine$double.eps^(1/12))
  expect_equal(as.numeric(s1), as.numeric(s2), tolerance = 1e-15)
  expect_error(suppressWarnings(Grad(x = 1:4, func = sum, method = "complex")),
               "Complex derivatives not implemented")
  expect_equal(as.numeric(Grad(x = 1:4, FUN = sum, side = NULL)), rep(1, 4), tolerance = 1e-10)
  # TODO: Richardson
})

test_that("gradients are correct", {
  f <- function(x) sum(sin(x))
  expect_equal(as.numeric(Grad(x = 1:4, f)), cos(1:4), tolerance = 1e-10)

  x <- structure(1:3, names = c("A", "B", "C"))
  g <- Grad(f, x, h = 0.01)
  expect_equal(attr(g, "step.size.method"), "user-supplied")
  expect_equal(names(g), c("A", "B", "C"))
})

test_that("Gradient step is auto-selected well", {
  f <- function(x) sum(sin(x))
  expect_equal(attr(Grad(x = 1:3, f), "step.size.method"), "default")
  expect_equal(attr(Grad(x = 1:3, f, h = 0.01), "step.size.method"), "user-supplied")
  expect_equal(attr(Grad(x = 1:3, f, h = "CR"), "step.size.method"), "CR")
  expect_equal(attr(Grad(x = 1:3, f, h = "SW"), "step.size.method"), "SW")
  g <- Grad(x = 1:3, FUN = f, h = "SW", elementwise = FALSE,
            vectorised = FALSE, multivalued = FALSE)
  expect_equal(attr(g, "step.search")[["exitcode"]], rep(0, 3), tolerance = 1e-15)
  expect_length(attr(g, "step.search")[["iterations"]], 3)
})

test_that("parallelisation of Grad works", {
  expect_equal(Grad(x = 1:3, FUN = sum, cores = 1),
               Grad(x = 1:3, FUN = sum, cores = 2))
})

test_that("function dimension check works", {
  f <- function(x) c(sum(sin(x)), sum(exp(x)))
  expect_error(Grad(f, 1:3), "for vector-valued functions")
  expect_error(Grad(f, 1:3, h = "SW"), "for vector-valued functions")
})

test_that("Grad works with automatic step sizes", {
  f <- function(x) x^2 - 2*x + 2
  expect_equal(as.numeric(Grad(f, x = 0.75, h = "CRm")), -0.5, tolerance = 1e-8)
})

test_that("Grad can accept dot arguments", {
  # These dot arguments must be accepted and used by checkDimensions()
  # and gradstep()
  f <- function(x, a0) sin(x + a0)
  expect_equal(Grad(f, x = 0, a0 = 1, h = 1e-5), Grad(sin, x = 1, h = 1e-5))
  expect_equal(attr(Grad(f, x = 0, a0 = 1, h = "SW"), "step.search")$exitcode, 0)
})

test_that("Grad can work on an arbitrary stencil", {
  d2 <- Grad(sin, 1, h = 1e-4, stencil = c(-1, 1))
  d4 <- Grad(sin, 1, h = 1e-4, stencil = c(-2, -1, 1, 2))
  expect_lt(abs(cos(1) - d4), abs(cos(1) - d2))
})
