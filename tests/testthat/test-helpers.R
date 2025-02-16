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

