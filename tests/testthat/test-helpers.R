test_that("Cluster creation handles inputs correctly", {
  expect_error(checkOrCreateCluster("rubbish"), "must be either")
  expect_error(checkOrCreateCluster(1:3), "is not a cluster")

  expect_equal(checkOrCreateCluster("lapply"), "lapply")
  expect_equal(checkOrCreateCluster("mclapply 1"), "lapply")

  cl <- parallel::makePSOCKcluster(2)
  expect_true(inherits(checkOrCreateCluster(cl), "cluster"))
  parallel::stopCluster(cl)
})

test_that("Parallel environments are set up correctly", {
  cl <- parallel::makePSOCKcluster(2)
  x <- setupParallelEnv(matrix, cl = cl, nrow = 3, data = 1:9)
  expect_true(is.environment(x$e))
  expect_equal(names(x$e), c("FUN", "data", "nrow"))
  expect_equal(x$FUN, x$e$FUN)
  parallel::stopCluster(cl)
})

test_that("Parallel runs are executed correctly", {
  cl <- parallel::makePSOCKcluster(2)
  x <- runParallel(sin, 1:3, cl = cl)
  expect_equal(x, as.list(sin(1:3)))
  expect_error(runParallel(sin, list(expression("sin(1)"))), "lapply should call the function directly")
  expect_error(runParallel(sin, list(expression("sin(1)")), cl = "mclapply 2"), "mclapply should call the function directly")
  expect_error(runParallel(sin, list(expression("sin(1)")), cl = "rubbish"), "should be")
  parallel::stopCluster(cl)
})

test_that("getExpr returns a valid expression", {
  e <- getExpr(327, rnorm, dots = list(sd = 0.1), fname = "rnorm")
  expect_true(grepl("rnorm\\(z, *sd *= *sd\\)", as.character(e)))
})

