test_that("input validation", {
  expect_error(Grad(x = 1:4, FUN = "rubbish"), "must be a function")
  expect_error(Grad(x = 1:4, FUN = sin, h = c(0.01, 0.02)), "must have length")
  expect_error(Grad(x = 1:4, FUN = sin, side = c(0, 1, 2, -2)), "must be 0 for central")
  expect_error(Grad(x = 1:4, FUN = sin, h = 0), "must be positive")
  expect_error(Grad(x = 1:4, FUN = sin, h = -0.001), "must be positive")
  w <- capture_warnings(Grad(x = 1:4, FUN = function(x) "0.1"))
  expect_true(any(grepl("at least one finite numeric", w)))
  expect_true(any(grepl("must output numeric values only", w)))
  expect_error(Grad(x = 1:4, FUN = sin, side = c(-1, 1)), "'side' argument must")
  expect_error(Grad(x = 1:4), "Pass the function")
  expect_error(Grad(x = 1:4, FUN = sin, deriv.order = c(1, 2)), "'deriv.order' must have length")
  expect_error(Grad(x = 1:4, FUN = sin, acc.order = c(1, 2)), "'acc.order' must have length")
  expect_warning(Grad(1:4, sin), "argument order")
})

test_that("vectorisation in GenD works", {
  expect_length(Grad(sin, 1:4), 4)
})

test_that("compatibility with numDeriv works", {
  expect_warning(Grad(sin, 1:4, method = "Richardson"), "numDeriv-like syntax")
})

test_that("missing values are treated properly", {
  noAttr <- function(x) { # Removing everything except for names
    attributes(x) <- attributes(x)["names"]
    x
  }
  f <- function(x) ifelse(abs(x-2) < 1e-3, x^2, NA)
  expect_equal(suppressWarnings(unname(noAttr(Grad(f, 1:3, side = 1)))), c(NA, 4, NA))
  expect_warning(Grad(f, 1:3, side = 1), "some non-numeric")
  f <- function(x) if (abs(x-2) < 1e-3) x^2 else "falsy"
  expect_warning(Grad(f, 1:3, side = 1, elementwise = TRUE, vectorised = FALSE,
                      multivalued = FALSE), "not even NA")
})
