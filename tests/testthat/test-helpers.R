test_that("Cores are checked well", {
  checkCores()
  checkCores(2)
  expect_warning(checkCores(1000), "more cores than")
})


test_that("Parallel runs are executed correctly", {
  cl <- parallel::makePSOCKcluster(2)
  x <- runParallel(FUN = sin, x = 1:3, cl = cl)
  expect_identical(x, as.list(sin(1:3)))
  parallel::stopCluster(cl)

  expect_error(runParallel(sin, 1:3, cl = "rubbish"), "passed as a cluster")
  expect_error(runParallel(sin, 1:3, cores = -1), "or a cluster object")
})


test_that("step size is valid for various values of x", {
  h0 <- stepx(10^(-16:3))
  h <- attr(GenD(sin, x = 10^(-16:3)), "step.size")
  expect_identical(h0, h)
  expect_true(all(diff(h) >= 0))
  expect_identical(h[1], h[2])  # Constant step for small x
  expect_identical(h[length(h)], 1000*.Machine$double.eps^(1/3))

  expect_true(all(diff(stepx(10^(-16:3), deriv.order = 2, acc.order = 4)) >= 0))
})


test_that("duplicated row indices are correct", {
  expect_identical(dupRowInds(mtcars[rep(1:10, 10), rep(1:10, 10)]), rep(1:10, 10))
})

test_that("Strings are aligned correctly", {
  x <- structure(1:4, names = month.name[1:4])
  o <- alignStrings(x, pad = "c")
  expect_identical(apply(nchar(o), 2, diff), rep(0L, 4))

  x <- matrix(c(1, 2.3, 4.567, 8, 9, 0), nrow = 2, byrow = TRUE)
  colnames(x) <- c("Andy", "Bradley", "Ci")
  o <- alignStrings(x, pad = "r")
  expect_identical(apply(nchar(o), 2, diff), matrix(0L, nrow = 2, ncol = 3))
})

test_that("Matrices are formatted and printed correctly", {
  x <- matrix(c(1234567, 12345.67, 123.4567, 1.23456, -1.23456e-1, 0,
                -1.23456e-4, 1.23456e-2, -1.23456e-6), nrow = 3)
  xf <- matrix(c("1.2e+6", " 1.2", "-1.2e-4", "1.2e+4", "-0.1", " 1.2e-2",
                 "1.2e+2", " 0  ", "-1.2e-6"), nrow = 3, byrow = TRUE)
  expect_identical(formatMat(x, digits = 1), xf)

  o <- printMat(x, digits = 1, shave.spaces = TRUE, begin = "c(", sep = ", ",
                end = ")", print = FALSE, format = TRUE)
  expect_identical(o[3], "c(1.2e+2,  0  , -1.2e-6)")
})
