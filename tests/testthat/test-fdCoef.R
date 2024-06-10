test_that("deriv.order >= 0 and acc.order > 0 must be integers", {
  expect_error(fdCoef(deriv.order = -1), "non-negative")
  expect_error(fdCoef(acc.order = 0), "positive")
})

test_that("stencil shorter than the needed order", {
  expect_error(fdCoef(2, stencil = 1), "stencil points")
  expect_warning(fdCoef(acc.order = 2, stencil = c(-2, 1)),
                 "achieve the requested accuracy")
})

test_that("side must be -1, 0, or 1", {
  expect_error(fdCoef(side = 3), "-1, 0, or 1")
  expect_warning(fdCoef(side = 2), "two-sided")
})

test_that("zero.action handling", {
  s <- -5:5
  expect_error(fdCoef(zero.action = "omit"), "zero.action")
  expect_length(fdCoef(stencil = s, zero.action = "none")$weights, 11)
  expect_length(fdCoef(stencil = s, zero.action = "drop")$weights, 10)
  expect_length(fdCoef(stencil = s, zero.action = "round")$weights, 11)
  expect_identical(unname(fdCoef(stencil = s, zero.action = "round")$weights[6]), 0)
})

test_that("correct weights for several stencils", {
  s05 <- c(-1, -0.5, 0.5, 1)
  expect_equal(fdCoef()$stencil, c(-1, 1))
  expect_equal(fdCoef()$weights, c(`x-1h` = -0.5, `x+1h` = 0.5))
  expect_equal(fdCoef(2), fdCoef(deriv.order = 2))
  expect_equal(unname(fdCoef(acc.order = 4)$weights) * 12, c(1, -8, 8, -1))
  expect_equal(unname(fdCoef(stencil = s05)$weights) * 6, c(1, -8, 8, -1))
  expect_equal(unname(fdCoef(3, stencil = s05)$weights), c(-4, 8, -8, 4))
  expect_equal(unname(fdCoef(stencil = 1:4)$weights), c(-13/3, 19/2, -7, 11/6))
  expect_warning(fdCoef(side = 0, acc.order = 1), "minimal stencil")
  expect_equal(attr(suppressWarnings(fdCoef(side = 0, acc.order = 1)), "accuracy.order"),
               c(requested = 1, effective = 2))
  expect_equal(attr(fdCoef(deriv.order = 3, stencil = -4:4), "accuracy.order"),
               c(requested = NA, effective = 6))
  expect_warning(fdCoef(stencil = c(-1, -1, 0, 1)), "duplicates")
  expect_equal(attr(suppressWarnings(fdCoef(stencil = c(-1, -1, 0, 1))),
                    "accuracy.order"), c(requested = NA, effective = 2))
  expect_warning(fdCoef(stencil = c(-0.002, -0.001, 0.001, 0.002)), "very close")
})
