test_that("input validation", {
  expect_error(Grad(x = 1:4, FUN = "rubbish"), "must be a function")
  expect_error(Grad(x = 1:4, FUN = sin, h = c(0.01, 0.02)), "must have length")
  expect_error(Grad(x = 1:4, FUN = sin, side = c(0, 1, 2, -2)), "must be 0 for central")
  expect_error(Grad(x = 1:4, FUN = sin, h = 0), "must be positive")
  expect_error(Grad(x = 1:4, FUN = sin, h = -0.001), "must be positive")
  expect_error(Grad(x = 1:4, FUN = function(x) "0.1"), "at least one finite numeric")
  expect_error(Grad(x = 1:4, FUN = sin, side = c(-1, 1)), "'side' argument must")
  expect_error(Grad(x = 1:4), "Pass the function")
  expect_error(Grad(x = 1:4, FUN = sin, deriv.order = c(1, 2)), "'deriv.order' must have length")
  expect_error(Grad(x = 1:4, FUN = sin, acc.order = c(1, 2)), "'acc.order' must have length")
  expect_warning(Grad(1:4, sin), "argument order")
})

test_that("vectorisation in GenD works", {
  expect_length(Grad(sin, 1:4), 4)
})
