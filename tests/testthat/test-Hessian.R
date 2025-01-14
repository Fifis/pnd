test_that("Hessians are correct", {
  x <- 1:4
  f <- function(x) prod(sin(x))
  fij <- Vectorize(function(i, j) prod(sin(x[-c(i, j)])) * prod(cos(x[c(i, j)])))
  h <- function(x) {
    m <- outer(seq_along(x), seq_along(x), fij)
    diag(m) <- -prod(sin(x))
    m
  }
  hes <- Hessian(f, x, report = 0)
  true.hes <- h(x)
  expect_equal(isSymmetric(hes), TRUE)
  expect_equal(hes, true.hes, tolerance = 1e-7)

  expect_equal(attr(Hessian(f, x), "step.size.method"), "default")
  expect_equal(attr(Hessian(f, x, h = 0.01), "step.size.method"), "user-supplied")
})

test_that("named arguments are handled correctly", {
  x <- 1:4
  f <- function(x) prod(sin(x))
  fij <- Vectorize(function(i, j) prod(sin(x[-c(i, j)])) * prod(cos(x[c(i, j)])))
  h <- function(x) {
    m <- outer(seq_along(x), seq_along(x), fij)
    diag(m) <- -prod(sin(x))
    m
  }
  hes <- Hessian(f, x, report = 0)
  true.hes <- h(x)

  names(x) <- LETTERS[1:4]
  hes.named <- Hessian(f, x, report = 0)
  expect_equal(colnames(hes.named), names(x))
  expect_equal(rownames(hes.named), names(x))

  x2 <- x
  names(x2) <- LETTERS[1:4]
  f2 <- function(x) prod(sin(x[c("A", "B")]) * sin(x[c("C", "D")]))
  hes2 <- Hessian(f2, x2, report = 0)
  expect_equal(unname(hes2), true.hes, tolerance = 1e-7)
})

test_that("Higher-order accuracy works", {
  x <- 1:4
  f <- function(x) prod(sin(x))
  fij <- Vectorize(function(i, j) prod(sin(x[-c(i, j)])) * prod(cos(x[c(i, j)])))
  h <- function(x) {
    m <- outer(seq_along(x), seq_along(x), fij)
    diag(m) <- -prod(sin(x))
    m
  }

  true.hes <- h(x)
  hes2 <- Hessian(f, x)
  hes4 <- Hessian(f, x, acc.order = 4)
  hes6 <- Hessian(f, x, acc.order = 6)
  expect_lt(max(abs(true.hes - hes4)), max(abs(true.hes - hes2)))
  expect_lt(max(abs(true.hes - hes6)), max(abs(true.hes - hes4)))

  # The auto-selected step must increase
  expect_equal(all(attr(hes2, "step.size") < attr(hes4, "step.size")), TRUE)
  expect_equal(all(attr(hes4, "step.size") < attr(hes6, "step.size")), TRUE)
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
  expect_error(Hessian(1:4, f), "argument order")
})

test_that("compatibility with numDeriv", {
  expect_warning(Hessian(x = 1:4, func = sum), "Use the argument")
})

# TODO: parallelisation
