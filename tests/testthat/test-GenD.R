test_that("input validation", {
  expect_error(Grad(x = 1:4, FUN = "rubbish", vectorised = FALSE), "must be a function")
  expect_error(Grad(x = 1:4, FUN = sin, h = c(0.01, 0.02), vectorised = FALSE), "must have length")
  expect_error(Grad(x = 1:4, FUN = sin, side = c(0, 1, 2, -2), vectorised = FALSE), "must be 0 for central")
  expect_error(Grad(x = 1:4, FUN = sin, h = 0, vectorised = FALSE), "must be positive")
  expect_error(Grad(x = 1:4, FUN = sin, h = -0.001, vectorised = FALSE), "must be positive")
  expect_warning(Grad(x = 1:4, FUN = function(x) "0.1", vectorised = FALSE), "numeric values only")
  expect_error(Grad(x = 1:4, FUN = sin, side = c(-1, 1), vectorised = FALSE), "'side' argument must")
  expect_error(Grad(x = 1:4), "Pass the function")
  expect_error(Grad(x = 1:4, FUN = sin, deriv.order = c(1, 2), vectorised = FALSE), "'deriv.order' must have length")
  expect_error(Grad(x = 1:4, FUN = sin, acc.order = c(1, 2), vectorised = FALSE), "'acc.order' must have length")
  expect_warning(Grad(1:4, sin), "argument order")
})

test_that("vectorisation in GenD works", {
  expect_length(Grad(sin, 1:4, vectorised = TRUE), 4)
})
