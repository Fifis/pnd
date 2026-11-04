test_that("input validation for GenD", {
  expect_warning(GenD(1:4, sin), "argument order")
})

test_that("input validation", {
  expect_error(GenD(x = 1:4, FUN = "rubbish"), "must be a function")
  expect_error(GenD(x = 1:4, FUN = sin, h = c(0.01, 0.02)), "must have length")
  expect_error(GenD(x = 1:4, FUN = sin, side = c(0, 1, 2, -2)), "must be 0 for central")
  expect_error(GenD(x = 1:4, FUN = sin, h = 0), "must be positive")
  expect_error(GenD(x = 1:4, FUN = sin, h = -0.001), "must be positive")

  w <- capture_warnings(GenD(x = 1:4, FUN = function(x) "0.1"))
  expect_true(any(grepl("at least one finite numeric", w)))
  expect_true(any(grepl("must output numeric values only", w)))

  expect_error(GenD(x = 1:4, FUN = sin, side = c(-1, 1)), "'side' argument must")
  expect_error(GenD(x = 1:4), "Pass the function")
  expect_error(GenD(x = 1:4, FUN = sin, deriv.order = c(1, 2)), "'deriv.order' must have length")
  expect_error(GenD(x = 1:4, FUN = sin, acc.order = c(1, 2)), "'acc.order' must have length")
  expect_warning(GenD(1:4, sin), "argument order")
})

test_that("vectorisation in GenD works", {
  expect_length(GenD(sin, 1:4), 4)
})

test_that("compatibility with numDeriv works", {
  expect_warning(GenD(sin, 1:4, method = "Richardson"), "numDeriv-like syntax")
})

test_that("missing values are treated properly", {
  noAttr <- function(x) { # Removing everything except for names
    attributes(x) <- attributes(x)["names"]
    x
  }
  f <- function(x) ifelse(abs(x-2) < 1e-3, x^2, NA)
  expect_equal(suppressWarnings(unname(noAttr(GenD(f, 1:3, side = 1)))), c(NA, 4, NA), tolerance = 1e-8)
  expect_warning(Grad(f, 1:3, side = 1), "some non-numeric")
  f <- function(x) if (abs(x-2) < 1e-3) x^2 else "falsy"
  expect_warning(GenD(f, 1:3, side = 1, elementwise = TRUE, vectorised = FALSE,
                      multivalued = FALSE), "not even NA")
})

test_that("custom stencils respect each coordinate's derivative order", {
  x <- c(A = 1, B = 2)
  h <- c(1/32, 1/16)
  orders <- c(1, 2)

  d <- GenD(function(x) sum(x^2), x, deriv.order = orders, stencil = -1:1, h = h)
  expect_equal(as.numeric(d), c(2, 2), tolerance = 1e-10)
  expect_named(d, names(x))

  d <- GenD(function(x) x^4, x, deriv.order = orders, stencil = -2:2, h = h)
  expect_equal(as.numeric(d), c(4, 48), tolerance = 1e-10)

  f <- function(x) c(sum(x^2), sum(x^3), sum(x^4))
  d <- GenD(f, x, deriv.order = orders, stencil = -2:2, h = h)
  expect_equal(as.numeric(d), c(2, 3, 4, 2, 12, 48), tolerance = 1e-10)

  expect_error(GenD(function(x) sum(x^2), x, deriv.order = orders,
                     stencil = c(-1, 1), h = h), "unique stencil points")
})

test_that("element-wise batching preserves derivatives and attributes", {
  named.input <- FALSE
  calls <- 0L
  f <- function(z, scale0) {
    calls <<- calls + 1L
    named.input <<- named.input || !is.null(names(z))
    scale0*z^4
  }
  for (x in list(c(A = 1), c(A = 1, B = 2, C = 3))) {
    n <- length(x)
    args <- list(FUN = f, x = x, h = c(1/32, 1/16, 1/8)[seq_len(n)],
                 deriv.order = seq_len(n), side = c(-1, 0, 1)[seq_len(n)],
                 acc.order = c(2, 4, 2)[seq_len(n)], scale0 = 2,
                 elementwise = TRUE, multivalued = FALSE, cores = 1)
    for (stencil in list(NULL, -3:3)) {
      calls <- 0L
      batched <- do.call(GenD, c(args, list(stencil = stencil, vectorised = TRUE)))
      expect_identical(calls, 1L)
      scalar <- do.call(GenD, c(args, list(stencil = stencil, vectorised = FALSE)))
      expect_identical(batched, scalar)
    }
  }
  expect_false(named.input)
})

test_that("unsuitable batches fall back to scalar evaluations", {
  run <- function(f, vectorised) GenD(f, c(-1, 1), h = 1/32, elementwise = TRUE,
                                      vectorised = vectorised, multivalued = FALSE, cores = 1)
  functions <- list(
    function(z) {if (length(z) > 1) stop("batch failure"); z^2},
    function(z) {if (length(z) > 1) warning("batch warning"); z^2},
    function(z) if (length(z) > 1) 0 else z^2,
    function(z) if (length(z) > 1) rep("bad", length(z)) else z^2,
    function(z) if (length(z) > 1) matrix(z^2, nrow = 1) else z^2,
    function(z) if (length(z) > 1) rep(Inf, length(z)) else z^2)
  for (f in functions) expect_identical(run(f, TRUE), run(f, FALSE))

  f <- function(z) {warning("scalar warning"); z^2}
  expect_identical(capture_warnings(batched <- run(f, TRUE)),
                   capture_warnings(scalar <- run(f, FALSE)))
  expect_identical(batched, scalar)

  f <- function(z) if (any(z < 0)) rep(NA_real_, length(z)) else z^2
  expect_warning(batched <- run(f, TRUE), "NA or NaN")
  expect_warning(scalar <- run(f, FALSE), "NA or NaN")
  expect_identical(batched, scalar)
  expect_equal(as.numeric(batched), c(NA, 2), tolerance = 1e-8)
})

test_that("element-wise derivatives retain the parallel path", {
  args <- list(FUN = sin, x = c(A = 1, B = 2), h = 1/32,
               elementwise = TRUE, vectorised = TRUE, multivalued = FALSE)
  expect_identical(do.call(GenD, c(args, list(cores = 1))),
                   suppressWarnings(do.call(GenD, c(args, list(cores = 2)))))
})

test_that("independent coordinate derivatives are equal to their counterparts with aggregation", {
  x <- c(A = 1, B = 2, C = 3)
  h <- setNames(c(1/32, 1/16, 1/8), names(x))
  f <- function(z) c(fourth = sum(z^4), fifth = sum(z^5))
  settings <- list(list(), list(deriv.order = c(1, 2, 3), stencil = -2:2),
                   list(acc.order = c(2, 4, 6)))
  for (setting in settings) {
    args <- c(list(FUN = f, x = x, h = h, elementwise = FALSE,
                   vectorised = FALSE, multivalued = TRUE), setting)
    d <- do.call(GenD, args)
    expect_identical(dimnames(d), list(c("fourth", "fifth"), names(x)))
    expect_identical(attr(d, "step.size"), h)
    expect_identical(attr(d, "step.size.method"), "user-supplied")
    for (i in seq_along(x)) {
      one <- args
      one$FUN <- function(z) {y <- x; y[i] <- z; f(y)}
      one$x <- x[i]
      one$h <- h[i]
      for (a in c("deriv.order", "acc.order"))
        if (length(one[[a]]) > 1) one[[a]] <- one[[a]][i]
      expect_identical(as.vector(d[, i]), as.vector(do.call(GenD, one)))
    }
  }
})
