test_that("Plug-in step handles inputs well", {
  expect_error(step.plugin(sin, 1, range = c(0, 1)), "must be a positive vector of length 2")
})

test_that("Plug-in step checks for small truncation errors", {
  expect_equal(step.plugin(sin, pi/3)$exitcode, 0)
  expect_equal(step.plugin(sin, pi/4)$exitcode, 1)
  expect_equal(step.plugin(sin, pi/2)$exitcode, 1)
})

test_that("Plug-in step checks for resprecting the ranges", {
  expect_equal(step.plugin(sin, pi/3, range = c(1e-4, 1e-3))$exitcode, 3)
  expect_equal(step.plugin(sin, pi/3, range = c(1e-7, 1e-6))$exitcode, 3)
})

test_that("Parallelisation in plug-in step selection works", {
  expect_equal(step.plugin(sin, 1, cores = 1), step.plugin(sin, 1, cores = 2))
})

test_that("User request in plug-in step selection for fewer cores is honoured", {
  expect_equal(step.plugin(sin, 1, cores = 1), step.plugin(sin, 1, cores = 2))
  f <- function(x) {
    if (length(x) > 1) stop("No vectorisation")
    Sys.sleep(0.05)
    sin(x)
  }
  t1 <- system.time(step.plugin(f, 1, cores = 1))
  t2 <- system.time(step.plugin(f, 1, cores = 2))
  expect_lt(t2[3], t1[3])
  clus <- parallel::makePSOCKcluster(2)
  expect_equal(step.plugin(sin, 1, cores = 1), step.plugin(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
