test_that("input validation", {
  expect_error(Grad(x = 1:4, FUN = "rubbish"), "must be a function")
  expect_error(Grad(x = 1:4, FUN = sum, h = c(0.01, 0.02)), "must have length")
  expect_error(Grad(x = 1:4, FUN = sum, side = c(0, 1, 2, -2)), "must contain values")
  expect_error(Grad(x = 1:4, FUN = sum, h = 0), "must be positive")
  expect_error(Grad(x = 1:4, FUN = sum, h = -0.001), "must be positive")
  expect_error(Grad(x = 1:4, FUN = function(x) "0.1"), "numeric values only")
})
