test_that("input validation works", {
  expect_error(fdCoef(deriv.order = -1), "non-negative")
  expect_error(fdCoef(acc.order = 0), "positive")
  expect_error(fdCoef(2, stencil = 1), "stencil points")
  expect_warning(fdCoef(acc.order = 2, stencil = c(-2, 1)),
                 "achieve the requested accuracy")
  expect_error(fdCoef(side = 3), "-1, 0, or 1")
  expect_warning(fdCoef(side = 2), "two-sided")
  expect_error(fdCoef(deriv.order = NULL), "'deriv.order' argument")
  expect_error(fdCoef(deriv.order = c(1, 2)), "'deriv.order' argument")
  expect_error(fdCoef(acc.order = c(1, 2)), "'acc.order' argument")
  expect_error(fdCoef(acc.order = NULL), "'acc.order' argument")
})

test_that("zero.action handling", {
  s <- -5:5
  expect_error(fdCoef(zero.action = "omit"), "zero.action")
  expect_length(fdCoef(stencil = s, zero.action = "none")$weights, 11)
  expect_length(fdCoef(stencil = s, zero.action = "drop")$weights, 10)
  expect_length(fdCoef(stencil = s, zero.action = "round")$weights, 11)
  expect_identical(unname(fdCoef(stencil = s, zero.action = "round")$weights[6]), 0)
})

test_that("default minimal stencils", {
  s1 <- fdCoef(side = -1)
  expect_equal(s1$stencil, -2:0, tolerance = 1e-15)
  expect_equal(unname(s1$weights), c(0.5, -2, 1.5), tolerance = 1e-15)

  s2 <- fdCoef(side = 1, acc.order = 3)
  expect_equal(s2$stencil, 0:3, tolerance = 1e-15)
  expect_equal(unname(s2$weights) * 6, c(-11, 18, -9, 2), tolerance = 1e-15)
})


test_that("correct weights for several stencils", {
  s05 <- c(-1, -0.5, 0.5, 1)
  expect_identical(fdCoef()$stencil, c(-1L, 1L))
  expect_identical(fdCoef()$weights, c(`x-1h` = -0.5, `x+1h` = 0.5))
  expect_identical(fdCoef(2), fdCoef(deriv.order = 2))
  expect_equal(unname(fdCoef(acc.order = 4)$weights) * 12, c(1, -8, 8, -1), tolerance = 1e-12)
  expect_identical(unname(fdCoef(stencil = s05)$weights) * 6, c(1, -8, 8, -1), tolerance = 1e-12)
  expect_identical(unname(fdCoef(3, stencil = s05)$weights), c(-4, 8, -8, 4), tolerance = 1e-12)
  expect_identical(unname(fdCoef(stencil = 1:4)$weights), c(-13/3, 19/2, -7, 11/6), tolerance = 1e-12)

  s1 <- suppressWarnings(fdCoef(side = 0, acc.order = 1))
  expect_identical(attr(s1, "accuracy.order"), c(requested = 1, effective = 2))
  expect_warning(fdCoef(side = 0, acc.order = 1), "requested 1st")
  expect_warning(fdCoef(side = 0, acc.order = 5), "requested 5th")
  expect_identical(attr(fdCoef(deriv.order = 3, stencil = -4:4), "accuracy.order"),
               c(requested = NA, effective = 6L))

  expect_warning(fdCoef(stencil = c(-0.002, -0.001, 0.001, 0.002)), "very close")
})

test_that("stencils of inappropriate length are correctly handled", {
  expect_warning(fdCoef(stencil = c(-1, -1, 0, 1)), "duplicates")
  expect_identical(attr(suppressWarnings(fdCoef(stencil = c(-1, -1, 0, 1))),
                    "accuracy.order"), c(requested = NA, effective = 2L))
  expect_warning(fdCoef(acc.order = 8, stencil = c(-2, -1, 1, 2)), "only 4th")
})

test_that("interpolation (degree-0 derivative) works", {
  expect_identical(attr(fdCoef(0, stencil = c(0.1, 0.2, 0.4, 0.8, 0.9) - 2/3),
                    "accuracy.order"), c(requested = NA, effective = 5L))

  f0 <- fdCoef(0, stencil = c(0.1, 0.2, 0.4, 0.8, 0.9) - 0.1)
  expect_identical(f0$stencil, 0)
  expect_identical(f0$weights, c(x = 1))
  expect_identical(attr(f0, "expansion"), "f exactly")
  expect_identical(attr(f0, "accuracy.order"), c(requested = NA, effective = Inf))
  expect_identical(attr(f0, "remainder.coef"), 0)
})
