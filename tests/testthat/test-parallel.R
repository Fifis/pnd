test_that("parallel derivatives preserve numerical results", {
  previous.options <- options(pnd.parallel = "auto", pnd.timing = FALSE,
                              pnd.max.cores = 2L)
  previous.cluster <- parallel::getDefaultCluster()
  on.exit({
    parallel::setDefaultCluster(previous.cluster)
    options(previous.options)
    stopParallel()
  })
  parallel::setDefaultCluster(NULL)
  stopParallel()

  f <- function(x) sum(x^2)
  outputs <- function(x) c(square = sum(x^2), linear = sum(x))
  environment(f) <- environment(outputs) <- baseenv()
  serial <- Grad(f, 1:2, h = 1/32, cores = 1)

  expect_identical(Grad(f, 1:2, h = 1/32, cores = 2), serial)
  expect_identical(Grad(f, 1:2, h = 1/32, cores = 2, preschedule = FALSE), serial)
  expect_equal(as.numeric(Jacobian(outputs, 1:2, h = 1/32, cores = 2)), c(2, 1, 4, 1), tolerance = 1e-8)
})
