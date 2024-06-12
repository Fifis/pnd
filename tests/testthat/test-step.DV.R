test_that("Dumontet-Vignes step selection handles inputs well", {
  f <- function(x) return(NA)
  expect_error(step.DV(x = 2, f), "Could not compute the function value")
})

test_that("Dumontet-Vignes step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.DV(x = 2, f, diagnostics = TRUE)
  expect_equal(s$exitcode, 0)
  expect_lt(s$abs.error, 1e-6)
  expect_equal(s$value, 32, tolerance = 1e-8)
  u <- s$iterations$ratio[length(s$iterations$ratio)]
  u <- max(u, 1/u)
  # Stopping criterion
  expect_gte(u, 2)
  expect_lte(u, 15)
})

test_that("Tweaking the DV algorithm for noisiser functions", {
  f <- function(x) x^4
  s.perfect <- step.DV(x = 2, f, h0 = 1e-7, alpha = 1)
  s.noisy <- step.DV(x = 2, f, h0 = 1e-7, alpha = 2)
  expect_lt(s.perfect$par, s.noisy$par)
})

