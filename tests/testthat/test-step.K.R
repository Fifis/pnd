test_that("Kostyrka's step selection handles inputs well", {
  expect_error(step.K(sin, 1, range = c(0, 1), cores = 1), "must be a positive vector of length 2")
  expect_warning(step.K(sin, 1, range = c(1e-4, 1e-6), cores = 1), "was extended to")
})

test_that("Kostyrka's step selection behaves reasonably", {
  s <- step.K(x = pi/4, sin, cores = 1)
  expect_identical(s$exitcode, 0L)
  expect_equal(s$value, sqrt(2)/2, tolerance = 1e-10)
})

test_that("Kostyrka's truncation estimates reproduce the cubic error", {
  s <- step.K(function(z) z^3, 0)
  h <- s$iterations$h
  error <- s$iterations$est.error[, "trunc"]
  expect_equal(mean(error/h^2, na.rm = TRUE), 1, tolerance = 1e-10)
})

test_that("Kostyrka's higher-order estimates reproduce the polynomial error", {
  s <- step.K(function(z) z^7, 0, deriv.order = 3, acc.order = 4)
  h <- s$iterations$h
  error <- s$iterations$est.error[, "trunc"]
  # The six-node target has K=1/10 and f^(7)=5040
  expect_equal(mean(error/h^4, na.rm = TRUE), 504, tolerance = 1e-10)
})

test_that("Kostyrka's method retries only unresolved auxiliary estimates", {
  calls <- 0L
  f <- function(z) { calls <<- calls + 1L; z^5 - 5*z^3 }
  k <- step.K(f, 0, cores = 1)
  m <- step.M(function(z) z^5 - 5*z^3, 0, cores = 1)
  kt <- k$iterations$est.error[, "trunc"]
  mt <- m$iterations$est.error[, "trunc"]
  good <- is.finite(mt) & mt > 0
  expect_equal(kt[good], mt[good], tolerance = 1e-12)
  # At h=1 the accuracy-two auxiliary estimate cancels exactly; accuracy four resolves it
  expect_equal(c(mt[m$iterations$h == 1], kt[k$iterations$h == 1]), c(0, 5))
  expect_equal(calls, 2*k$counts)
})

test_that("Kostyrka interpolates the leading error at the returned step", {
  s <- step.K(function(z) z^3, 1)
  expect_equal(s$value, 3 + s$par^2, tolerance = 1e-12)
})

test_that("Kostyrka's fallback uses the model rounding target", {
  f <- function(z) ifelse(abs(z) %in% c(1/4, 1/2, 1), 1 + z^3, NA_real_)
  p <- 1e-8
  s <- suppressWarnings(step.K(f, 0, max.rel.error = p))
  # K*|f'''|=1 and the finest complete target stencil has maximum 65/64
  B <- (p + .Machine$double.eps/2) * 65/64
  expect_equal(s$iterations$fallback$target, (2*B^2)^(1/3), tolerance = 1e-12)
})

test_that("Kostyrka's step selection respects non-integer stencil offsets", {
  expect_equal(step.K(exp, 1, shrink.factor = 0.6, cores = 1)$value, exp(1), tolerance = 1e-7)
  expect_equal(step.K(exp, 1, shrink.factor = 0.8, cores = 1)$value, exp(1), tolerance = 1e-7)
  expect_equal(step.K(exp, 1, shrink.factor = 1/sqrt(2), cores = 1)$value, exp(1), tolerance = 1e-7)

  # Higher accuracy order also works
  expect_equal(suppressWarnings(step.K(identity, 1, acc.order = 4, shrink.factor = 0.6, cores = 1))$value, 1, tolerance = 1e-10)
})

test_that("Kostyrka's step selection works for gradients", {
  f <- function(x) sin(x[1]) * cos(x[2]) + exp(x[3])
  s <- gradstep(f, 1:3, method = "K")
  expect_identical(s$exitcode, rep(0L, 3))
})

test_that("Kostyrka's method returns reasonable values with unfortunate inputs", {
  s <- step.K(sin, 1, shrink.factor = 1/8, max.rel.error = 0.5, cores = 1)
  expect_equal(s$value, cos(1), tolerance = 1e-9)
  expect_identical(s$exitcode, 0L)

  f <- function(x) ifelse(x %in% 0:2, x^2, NA)
  expect_identical(step.K(f, 1, cores = 1)$exitcode, 2L)
})

test_that("Kostyrka's method returns reasonable values with tricky functions", {
  s <- step.K(function(x) x^3+1/x, 1)
  expect_identical(s$exitcode, 0L)
  # The leading third and fifth derivatives vanish, so a fallback is legitimate
  s <- suppressWarnings(step.K(function(x) 6*x^5 - 56*x^3 + 1/x^2, 1))
  expect_equal(s$value, -140, tolerance = 1e-10)
})

test_that("Kostyrka's algorithm returns expected non-zero exit codes", {
  # Very noisy right branch
  expect_identical(suppressWarnings(step.K(function(x) x^2, x = 1e-8)$exitcode), 1L)
  # No right branch at all
  s <- suppressWarnings(step.K(function(x) x^2, x = 0))
  expect_identical(s$exitcode, 2L)
  expect_null(s$iterations$fit)
})

test_that("Kostyrka's fallback uses a finite upper neighbour", {
  f <- function(z) ifelse(abs(z-1) < 1e-10 | abs(z-1) >= 0.1, z^3, NA_real_)
  s <- suppressWarnings(step.K(f, 1, range = c(2^-46, 2^14), cores = 1))
  # Only h = 1/8 has a finite derivative among the two neighbours
  expect_equal(s$value, 3 + (1/8)^2, tolerance = 1e-10)
})

test_that("Kostyrka's fallback uses a finite lower neighbour", {
  f <- function(z) {
    r <- abs(z-1)
    ifelse(r < 1e-10 | r >= 0.15 | (r >= 1/16 & r < 1/8), z^3, NA_real_)
  }
  s <- suppressWarnings(step.K(f, 1, range = c(2^-46, 2^14), cores = 1))
  # Only h = 1/16 has a finite derivative among the two neighbours.
  expect_equal(s$value, 3 + (1/16)^2, tolerance = 1e-10)
})

test_that("Kostyrka's fallback uses a derivative from the full grid", {
  f <- function(z) ifelse(abs(z-1) < 1e-10 | abs(z-1) >= 0.15, z^3, NA_real_)
  s <- suppressWarnings(step.K(f, 1, range = c(2^-46, 2^14), cores = 1))
  # Both neighbours are missing; h = 1/4 is the finite point nearest the grid's centre.
  expect_equal(s$value, 3 + (1/4)^2, tolerance = 1e-10)
})

test_that("Kostyrka's fallback fails explicitly without a finite derivative", {
  f <- function(z) ifelse(abs(z) < 1e-160, sign(z)*1e150, NA_real_)
  # Function values are finite on part of the grid, but every available derivative overflows
  expect_error(suppressWarnings(step.K(f, 0, range = c(1e-170, 1), cores = 1)), "finite derivative")
})

test_that("Parallelisation in Kostyrka's algorithm works", {
  expect_identical(step.K(sin, 1, cores = 1), step.K(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.K(sin, 1, cores = 1), step.K(sin, 1, cl = clus))
  parallel::stopCluster(clus)

  # Testing a slow function
  # f <- function(x) {Sys.sleep(0.1); sin(x)}
  # system.time(step.K(f, 1))
  # system.time(step.K(f, 1, cores = 12))
})
