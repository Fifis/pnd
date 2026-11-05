withDerivativeTiming <- function(code) {
  previous.options <- options(pnd.timing = TRUE)
  previous.cluster <- parallel::getDefaultCluster()
  on.exit({
    options(previous.options)
    parallel::setDefaultCluster(previous.cluster)
  })
  parallel::setDefaultCluster(NULL)
  force(code)
}

test_that("timing can be enabled and disabled", {
  withDerivativeTiming({
    timed <- Grad(sum, 1:2, h = 1/32)
    timing <- attr(timed, "timing")
    expect_named(timing, c("preparation", "evaluation", "elapsed"))
    expect_true(all(is.finite(timing) & timing >= 0))
    expect_equal(attr(timed, "function.calls"), 4)

    options(pnd.timing = FALSE)
    plain <- Grad(sum, 1:2, h = 1/32)
    expect_null(attr(plain, "timing"))
    expect_null(attr(plain, "function.calls"))
    expect_true(is.na(attr(checkDimensions(sum, 1:2), "seconds")))
    attr(timed, "timing") <- attr(timed, "function.calls") <- NULL
    expect_identical(timed, plain)
  })
})

test_that("Jacobian and Hessian return timing information", {
  withDerivativeTiming({
    jac <- Jacobian(function(x) c(sum(x), sum(x^2)), 1:2, h = 1/32)
    hes <- Hessian(function(x) sum(x^2), 1:2, h = 1/32)
    expect_length(attr(jac, "timing"), 3)
    expect_equal(attr(jac, "function.calls"), 4)
    expect_length(attr(hes, "timing"), 3)
    expect_equal(attr(hes, "function.calls"), 9)

    options(pnd.timing = FALSE)
    hes <- Hessian(sum, 1:2, h = 1/32)
    expect_null(attr(hes, "timing"))
    expect_null(attr(hes, "function.calls"))
  })
})

test_that("call counts include batch attempts and exclude preparation", {
  withDerivativeTiming({
    batch <- GenD(sin, 1:2, h = 1/32)
    expect_equal(attr(batch, "function.calls"), 1)
    f <- function(x) {if (length(x) != 1) stop("Scalar input required"); x^2}
    fallback <- GenD(f, 1:2, h = 1/32, elementwise = TRUE,
                     vectorised = TRUE, multivalued = FALSE)
    expect_equal(attr(fallback, "function.calls"), 5)

    calls <- 0L
    f <- function(x) {calls <<- calls + 1L; sum(x^2)}
    automatic <- Grad(f, 1:2, h = "SW")
    expect_gt(calls, 4)
    expect_equal(attr(automatic, "function.calls"), 4)
  })
})

test_that("PSOCK timings work with reused workers", {
  withDerivativeTiming({
    previous <- options(pnd.parallel = "PSOCK", pnd.max.cores = 2L)
    on.exit({options(previous); stopParallel()}, add = TRUE)
    stopParallel()

    Grad(sum, 1:2, h = 1/32, cores = 2)
    cl <- getParallelCluster(2)
    result <- Grad(sum, 1:2, h = 1/32, cores = 2)
    expect_equal(as.numeric(result), c(1, 1))
    expect_length(attr(result, "timing"), 3)
    expect_equal(attr(result, "function.calls"), 4)
    expect_identical(getParallelCluster(2), cl)

    stopParallel()
    expect_null(.pndCluster$cl)
  })
})
