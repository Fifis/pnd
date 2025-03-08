test_that("Cores are checked well", {
  checkCores()
  checkCores(2)
  expect_warning(checkCores(1000), "more cores than")
})


test_that("Parallel runs are executed correctly", {
  cl <- parallel::makePSOCKcluster(2)
  x <- runParallel(FUN = sin, x = 1:3, cl = cl)
  expect_equal(x, as.list(sin(1:3)))
  parallel::stopCluster(cl)

  expect_error(runParallel(sin, 1:3, cl = "rubbish"), "passed as a cluster")
  expect_error(runParallel(sin, 1:3, cores = -1), "or a cluster object")
})


test_that("step size is valid for various values of x", {
  h0 <- stepx(10^(-16:3))
  h <- attr(GenD(sin, x = 10^(-16:3)), "step.size")
  expect_equal(h0, h)
  expect_true(all(diff(h) >= 0))
  expect_equal(h[1], h[2])  # Constant step for small x
  expect_equal(h[length(h)], 1000*.Machine$double.eps^(1/3))

  expect_true(all(diff(stepx(10^(-16:3), deriv.order = 2, acc.order = 4)) >= 0))
})

