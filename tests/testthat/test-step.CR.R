test_that("Curtis-Reid step selection handles inputs well", {
  f <- function(x) return(NA)
  expect_error(step.CR(x = 2, f), "Could not compute the function value")
  expect_error(step.CR(sin, 1, version = "new"), "should be one of")
  expect_error(step.CR(sin, 1, tol = 1e-4), "must be >=")
  expect_warning(step.CR(sin, 1, acc.order = 4), "Using acc.order")
  expect_error(step.CR(sin, 1, range = c(0, 1)), "must be a positive vector of length 2")
})


test_that("Curtis-Reid step selection behaves reasonably", {
  f <- function(x) x^4
  s <- step.CR(x = 2, f)
  expect_identical(s$exitcode, 0L)
  expect_lt(sum(s$abs.error), 1e-4)
  expect_equal(s$value, 32, tolerance = 1e-8)
  u <- s$iterations$ratio[length(s$iterations$ratio)]
  expect_gt(u, 10)
  expect_lt(u, 1000)

  s2 <-  step.CR(x = 2, f, version = "modified")
  expect_lt(sum(s2$abs.error), 1e-7)
  expect_equal(s2$value, 32, tolerance = 1e-8)

  s3 <- step.CR(x = 2, f, version = "modified", acc.order = 4)
  expect_lt(sum(s2$abs.error), 5e-8)
  expect_equal(s3$value, 32, tolerance = 1e-8)

  s4 <- step.CR(x = sqrt(2), FUN = function(x) x^6 -2*x^4 - 4*x^2, h0 = 2^-16)
  expect_lt(abs(s4$value), 1e-8)

})

test_that("Curtis--Reid algorithm stops if the function returns NA for all allowed step sizes", {
  f <- function(x) ifelse(abs(x - 2) < 1e-8, x^4, NA)
  expect_error(step.CR(f, 2, range = c(1e-7, 1e-2)), "attempts of step shrinkage")
})

test_that("Curtis--Reid steps grow for linear functions", {
  expect_identical(step.CR(x = 1, function(x) x)$exitcode, 1L)
})

test_that("Large and small initial values in Curtis--Reid cause range problems", {
  f <- function(x) x^4
  expect_identical(step.CR(x = 2, f, h0 = 1000)$exitcode, 3L)
  expect_identical(step.CR(x = 2, f, h0 = 1000, range = c(1e-10, 1e5))$exitcode, 0L)
  expect_identical(step.CR(x = 2, f, h0 = 1e-12)$exitcode, 1L)  # In theory, it should be 3,
  # but the truncation error estimate is really zero...
  expect_identical(step.CR(x = 2, f, h0 = 9e-9, range = c(1e-8, 2e-8))$exitcode, 3L)
  expect_identical(step.CR(x = 2, FUN = f, h0 = 1000, maxit = 2, range = c(1e-12, 1))$exitcode, 4L)
})

# test_that("Parallelisation speeds things up for Curtis--Reid", {
#   f <- function(x) {Sys.sleep(0.3); return(sin(x))}
#   cl <- parallel::makeCluster(2)
#   t1 <- system.time(step.CR(f, pi/4, version = "modified", cores = 1))
#   t2 <- system.time(step.CR(f, pi/4, version = "modified", cl = cl))
#   parallel::stopCluster(cl)
#   expect_lt(t2[3], t1[3])
# })

test_that("Parallelisation in CR works", {
  expect_identical(step.CR(sin, 1, cores = 1), step.CR(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.CR(sin, 1, cores = 1), step.CR(sin, 1, cl = clus))
  parallel::stopCluster(clus)
})
