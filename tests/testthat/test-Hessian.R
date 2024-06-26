test_that("Hessians are correct", {
  x <- 1:4
  f <- function(x) prod(sin(x))
  h <- function(x) {
    n <- length(x)
    m <- outer(1:n, 1:n, Vectorize(function(i, j) {
      o <- sin(x)
      d <- cos(x)
      o[c(i, j)] <- d[c(i, j)]
      prod(o)
    }))
    diag(m) <- -prod(sin(x))
    m
  }
  hes <- Hessian(f, 1:4, report = 0)
  true.hes <- h(1:4)
  expect_equal(isSymmetric(hes), TRUE)
  expect_equal(hes, true.hes, tolerance = 1e-7)

  names(x) <- LETTERS[1:4]
  hes2 <- Hessian(f, x, report = 0)
  expect_equal(colnames(hes2), names(x))
  expect_equal(rownames(hes2), names(x))

  expect_equal(attr(Hessian(f, x), "step.size.method"), "default")
  expect_equal(attr(Hessian(f, x, h = 0.01), "step.size.method"), "user-supplied")
})

test_that("function dimension check works", {
  f <- function(x) c(sin(x), exp(x))
  expect_error(Hessian(f, 1:3), "only scalar functions")
})

test_that("input check works", {
  f <- function(x) prod(sin(x))
  expect_equal(Hessian(f, 1:4, side = NULL), Hessian(f, 1:4, side = 0))

  expect_error(Hessian(x = 1:3, FUN = "sin"), "must be a function")
  expect_error(Hessian(f, 1:4, side = 2), "'side' argument")
  expect_error(Hessian(as.character, 1:3), "numeric values only")
  expect_error(Hessian(f, 1:3, h = "SW"), "algorithms not implemented")
  expect_error(Hessian(f, 1:3, h = -1), "must be positive")
  expect_error(suppressWarnings(Hessian(as.character, 1:3, acc.order = 1:2)),
  "'acc.order' must have length")
  expect_error(suppressWarnings(Hessian(as.character, 1:3, h = 1:2)),
               "must have length")
})

test_that("compatibility with numDeriv", {
  expect_warning(Hessian(x = 1:4, func = sum), "Use the argument")
})

# TODO: parallelisation
