test_that("Correct weights for several stencils", {
  expect_equal(fdCoef(),
               list(stencil = c(-1, 1),
                    weights = c(`x-1h` = -0.5, `x+1h` = 0.5)),
               tolerance = 1e-15)
})
