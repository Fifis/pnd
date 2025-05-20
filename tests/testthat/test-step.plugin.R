test_that("Plug-in step handles inputs well", {
  expect_error(step.plugin(sin, 1, range = c(0, 1)), "must be a positive vector of length 2")
})

test_that("Plug-in step checks for small truncation errors", {
  expect_identical(step.plugin(sin, pi/3)$exitcode, 0L)
  expect_identical(step.plugin(function(x) x, pi/2)$exitcode, 1L)
})

test_that("Plug-in step checks for resprecting the ranges", {
  expect_identical(step.plugin(sin, pi/3, range = c(1e-4, 1e-3))$exitcode, 3L)
  expect_identical(step.plugin(sin, pi/3, range = c(1e-9, 1e-8))$exitcode, 3L)
})

test_that("Plug-in algorithm stops if the function returns NA for all allowed step sizes", {
  f <- function(x) ifelse(abs(x - 2) < 1e-8, x^4, NA)
  expect_error(step.plugin(f, 2, range = c(1e-7, 1e-2)), "attempts of step shrinkage")
})

test_that("Parallelisation in plug-in step selection works", {
  expect_identical(step.plugin(sin, 1, cores = 1), step.plugin(sin, 1, cores = 2))
})

test_that("User request in plug-in step selection for fewer cores is honoured", {
  expect_identical(step.plugin(sin, 1, cores = 1), step.plugin(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.plugin(sin, 1, cores = 1), step.plugin(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
