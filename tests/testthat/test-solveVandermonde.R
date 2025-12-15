test_that("Stencil contains unique values", {
  expect_error(solveVandermonde(s = c(1, 1, 2), b = 1:3),
               "must be unique")
})

test_that("Stencil and coefficients are the same length", {
  expect_error(solveVandermonde(s = 0:1, b = 1:3),
               "must be of the same length")
})

test_that("Unsorted stencil are bad", {
  expect_warning(solveVandermonde(s = c(0, 1, -2), b = c(0, 0, 1)),
                 "unsorted")
})

test_that("Vandermonde system for positive stencils", {
  s <- 0:5
  b <- (s == 4) * 24
  result <- solveVandermonde(s, b)
  expect_type(result, "double")
  expect_length(result, length(s))
})
