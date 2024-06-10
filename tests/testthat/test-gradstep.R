test_that("step size length h0 must be 1 or length(x)", {
  expect_error(gradstep(x = 1, FUN = sin, h0 = c(1e-5, 1e-6)), "must be a scalar")
  expect_error(gradstep(x = 1:2,
                      FUN = function(z) {if (length(z)>1) stop("Non-vectorised"); z^2}),
               "must be finite")
})

test_that("method arguments with strange names cause an error", {
  expect_error(gradstep(x = 1, FUN = sin, method = "CR", method.args = list(rubbish = TRUE)),
               "arguments are not supported")
})

test_that("some methods do not converge for unfortunate inputs", {
  expect_equal(gradstep(x = 1e10, FUN = sin, h0 = 1e-20, method = "SW")$exitcode, 4)
  expect_equal(gradstep(x = 1, FUN = sin, method = "DV",
                      method.args = list(range = c(1e-20, 1e-22)))$exitcode, 0)

})
