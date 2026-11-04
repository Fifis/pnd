withDerivativeTiming <- function(code, enabled = TRUE) {
  previous.options <- options(pnd.timing = enabled)
  previous.cluster <- parallel::getDefaultCluster()
  on.exit({
    options(previous.options)
    parallel::setDefaultCluster(previous.cluster)
  })
  parallel::setDefaultCluster(NULL)
  force(code)
}

withoutDerivativeTiming <- function(x) {
  for (a in c("timing", "function.calls")) attr(x, a) <- NULL
  x
}

expectDerivativeTiming <- function(x, calls) {
  timing <- attr(x, "timing")
  expect_named(timing, c("preparation", "evaluation", "elapsed"))
  expect_equal(attr(x, "function.calls"), calls)
}

test_that("enabled timing attaches the timing attribute", {
  previous <- options(pnd.timing = TRUE)
  on.exit(options(previous))
  result <- Grad(sum, 1:3)
  expect_true("timing" %in% names(attributes(result)))
})

test_that("timing can be switched off", {
  previous <- options(pnd.timing = TRUE)
  on.exit(options(previous))
  result <- Grad(sum, 1:3)
  expect_true("timing" %in% names(attributes(result)))
  options(pnd.timing = FALSE)
  result <- Grad(sum, 1:3)
  expect_null(attr(result, "timing"))
  expect_null(attr(result, "function.calls"))
  expect_true(is.na(attr(checkDimensions(sum, 1:3), "seconds")))
})

test_that("callbacks preserve lazy, quoted, and timing arguments", {
  withDerivativeTiming({
    f <- function(z, unused, expr, timing) {
      stopifnot(is.call(expr))
      timing * sum(z^2)
    }
    d <- GenD(f, 1:3, h = 1/32, elementwise = FALSE,
              vectorised = FALSE, multivalued = FALSE,
              unused = stop("Unused argument was evaluated"),
              expr = quote(stop("Quoted call was evaluated")), timing = 2)
    expect_equal(as.numeric(d), c(4, 8, 12))
    expectDerivativeTiming(d, 6)
  })
})

test_that("Jacobian timings honour supplied initial values", {
  withDerivativeTiming({
    f <- function(z) c(square = sum(z^2), linear = sum(z), cube = sum(z^3))
    automatic <- Jacobian(f, 1:3, h = 1/32)
    supplied <- Jacobian(f, 1:3, h = 1/32, f0 = f(1:3))
    expectDerivativeTiming(automatic, 6)
    expectDerivativeTiming(supplied, 6)
    expect_identical(withoutDerivativeTiming(automatic), withoutDerivativeTiming(supplied))
  })
})

test_that("batch trials and fallback evaluations are counted", {
  withDerivativeTiming({
    args <- list(x = 1:3, h = 1/32, elementwise = TRUE,
                 vectorised = TRUE, multivalued = FALSE)
    expectDerivativeTiming(do.call(GenD, c(list(FUN = sin), args)), 1)
    f <- function(z) {if (length(z) != 1L) stop("Scalar input required"); z^2}
    batched <- do.call(GenD, c(list(FUN = f), args))
    expectDerivativeTiming(batched, 7)
    args$vectorised <- FALSE
    scalar <- do.call(GenD, c(list(FUN = f), args))
    expectDerivativeTiming(scalar, 6)
    expect_identical(withoutDerivativeTiming(batched), withoutDerivativeTiming(scalar))
  })
})

test_that("step selection is preparation and final stencil calls are counted", {
  withDerivativeTiming({
    calls <- 0L
    f <- function(z) {calls <<- calls + 1L; sum(z^2)}
    d <- Grad(f, 1:3, h = "SW")
    expect_gt(calls, 6)
    expectDerivativeTiming(d, 6)
  })
})

test_that("timing preserves warnings and failed evaluation results", {
  withDerivativeTiming({
    args <- list(x = 1:3, h = 1/32, elementwise = TRUE,
                 vectorised = TRUE, multivalued = FALSE)
    f <- function(z) {warning("Evaluation warning"); z^2}
    warnings <- capture_warnings(timed <- do.call(GenD, c(list(FUN = f), args)))
    expectDerivativeTiming(timed, 7)
    options(pnd.timing = FALSE)
    plain.warnings <- capture_warnings(plain <- do.call(GenD, c(list(FUN = f), args)))
    expect_identical(warnings, plain.warnings)
    expect_length(warnings, 6)
    expect_identical(withoutDerivativeTiming(timed), plain)
    options(pnd.timing = TRUE)
    f <- function(z) stop("Evaluation failed")
    args$elementwise <- args$vectorised <- FALSE
    expect_warning(failed <- do.call(GenD, c(list(FUN = f), args)), "NA or NaN")
    expectDerivativeTiming(failed, 6)
    expect_true(all(is.na(failed)))
  })
})

test_that("nested derivatives return independent measurements", {
  withDerivativeTiming({
    inner <- NULL
    f <- function(z) {
      inner <<- Grad(function(y) sum(y^2), 1:2, h = 1/32)
      sum(z^2)
    }
    outer <- GenD(f, 1:3, h = 1/32, elementwise = FALSE,
                  vectorised = FALSE, multivalued = FALSE)
    expectDerivativeTiming(outer, 6)
    expectDerivativeTiming(inner, 4)
  })
})

test_that("forked evaluations return timings under both schedules", {
  skip_if(.Platform$OS.type == "windows", "Forking is unavailable on Windows")
  withDerivativeTiming({
    args <- list(FUN = sum, x = 1:3, h = 1/32)
    serial <- do.call(Grad, args)
    for (schedule in c(TRUE, FALSE)) {
      forked <- do.call(Grad, c(args, list(cores = 2, preschedule = schedule)))
      expectDerivativeTiming(forked, 6)
      expect_identical(withoutDerivativeTiming(forked), withoutDerivativeTiming(serial))
    }
  })
})

test_that("reused explicit and default PSOCK clusters preserve values", {
  withDerivativeTiming({
    cl <- parallel::makePSOCKcluster(2)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    setup <- function(paths) {.libPaths(paths); loadNamespace("pnd"); NULL}
    environment(setup) <- baseenv()
    parallel::clusterCall(cl, setup, .libPaths())
    args <- list(FUN = sum, x = 1:3, h = 1/32)
    serial <- do.call(Grad, args)
    jargs <- list(FUN = function(z) c(square = sum(z^2), cube = sum(z^3)),
                  x = 1:3, h = 1/32, elementwise = FALSE,
                  vectorised = FALSE, multivalued = TRUE)
    environment(jargs$FUN) <- baseenv()
    serial.j <- do.call(Jacobian, jargs)
    for (schedule in c(TRUE, FALSE)) {
      explicit <- do.call(Grad, c(args, list(cl = cl, preschedule = schedule)))
      explicit.j <- do.call(Jacobian, c(jargs, list(cl = cl, preschedule = schedule)))
      expectDerivativeTiming(explicit.j, 6)
      expect_identical(withoutDerivativeTiming(explicit.j), withoutDerivativeTiming(serial.j))
      parallel::setDefaultCluster(cl)
      default <- do.call(Grad, c(args, list(preschedule = schedule)))
      parallel::setDefaultCluster(NULL)
      for (d in list(explicit, default)) {
        expectDerivativeTiming(d, 6)
        expect_identical(withoutDerivativeTiming(d), withoutDerivativeTiming(serial))
      }
    }
    lazy <- function(z, unused) sum(z^2)
    environment(lazy) <- baseenv()
    d <- GenD(lazy, 1:3, h = 1/32, elementwise = FALSE,
              vectorised = FALSE, multivalued = FALSE, cl = cl,
              unused = stop("Unused cluster argument was evaluated"))
    expectDerivativeTiming(d, 6)
    expect_equal(as.numeric(d), c(2, 4, 6))
  })
})
