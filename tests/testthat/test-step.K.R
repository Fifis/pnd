test_that("Kostyrka's step selection handles inputs well", {
  expect_error(step.K(sin, 1, range = c(0, 1), cores = 1), "must be a positive vector of length 2")
  expect_warning(step.K(sin, 1, range = c(1e-4, 1e-6), cores = 1), "was extended to")
})

test_that("Kostyrka's step selection behaves reasonably", {
  s <- step.K(x = pi/4, sin, cores = 1)
  expect_identical(s$exitcode, 0L)
  expect_equal(s$value, sqrt(2)/2, tolerance = 1e-10)
})

test_that("Unresolved affine estimates cannot create a rounding branch", {
  calls <- 0L
  f <- function(z) { calls <<- calls + 1L; z }
  s <- step.K(f, 1, cores = 1)
  error <- s$iterations$est.error[, "trunc"]
  expect_equal(c(s$value, calls, sum(is.finite(error) & error > 0)), c(1, 2*s$counts, 0))
})

test_that("Changing function units preserves unresolved estimates", {
  base <- step.K(sin, 1, cores = 1)
  scaled <- step.K(function(z) 2^40*sin(z), 1, cores = 1)
  expect_identical(base$iterations$est.error[, "trunc"] == 0,
                   scaled$iterations$est.error[, "trunc"] == 0)
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

test_that("Kostyrka's method prefers the wider auxiliary stencil", {
  k <- step.K(function(z) z^5 - 5*z^3, 0, cores = 1)
  # The third derivative at zero is -30; the shorter formula cancels at h=1
  expect_equal(unname(k$iterations$est.error[k$iterations$h == 1, "trunc"]), 5)
  expect_identical(k$iterations$aprime[k$iterations$h == 1], 4L)
})

test_that("Kostyrka's algorithm retries unresolved wider estimates without new calls", {
  calls <- 0L
  # At h=1 the wider third-derivative formula cancels for this polynomial
  f <- function(z) { calls <<- calls + 1L; z^7 + 84*z^3 }
  k <- step.K(f, 0, cores = 1)
  expect_identical(k$iterations$aprime[k$iterations$h == 1], 2L)
  expect_gt(k$iterations$est.error[k$iterations$h == 1, "trunc"], 0)
  expect_equal(calls, 2*k$counts)
})

test_that("Kostyrka's method reduces cancellation when the truncation branch vanishes", {
  s <- suppressWarnings(step.K(function(z) pi*z + 2, 1, max.rel.error = .Machine$double.eps/2))
  expect_lt(abs(s$value/pi - 1), 16*.Machine$double.eps)
  expect_equal(s$par, s$iterations$fallback$stable.step)
})

test_that("Kostyrka checks finer estimates before trusting a coarse plateau", {
  f <- function(z) z + (z-1)*(1-exp(-((z-1)/.003)^4))
  s <- suppressWarnings(step.K(f, 1))
  # The derivative is 1 at the centre and approaches 2 outside the narrow transition
  expect_equal(s$value, 1, tolerance = 1e-7)
})

test_that("Kostyrka checks higher-order error when the third derivative vanishes", {
  f <- function(z) 1e6 + .125*(z-12) + .001*(z-12)^5
  s <- suppressWarnings(step.K(f, 12))
  # A zero third derivative still leaves an h^4 error in the target formula
  expect_lt(abs(s$value - .125), 1e-7)
})

test_that("Kostyrka interpolates the leading error at the returned step", {
  s <- step.K(function(z) z^3, 0)
  expect_equal(s$value, s$par^2, tolerance = 1e-12)
})

test_that("Kostyrka balances truncation and interpolation noise", {
  # D(1) and D(2) have bias 1 and 4, and independent noise variances 64 and 16
  # Minimising their mixture's squared error gives weight 61/89 on D(2)
  noise <- c(ratio = 1, overlap = 0)
  step <- interpolationStepK(c(1, 2), c(TRUE, TRUE), 2, 1, 1, 2, noise)
  expect_equal(step^2, 272/89)
  expect_equal(interpolationStepK(c(1, 2), c(TRUE, TRUE), 1.25, 1, 1, 2, noise), 1.25)
})

test_that("Kostyrka accounts for the shared centre of even derivatives", {
  # Second differences share their centre: its squared weight is 4 of the total 6
  # Biases 1 and 4, variances 16 and 1, and covariance 8/3 give equal optimal weights
  noise <- c(ratio = 1, overlap = 2/3)
  step <- interpolationStepK(c(1, 2), c(TRUE, TRUE), 2, .5, 2, 2, noise)
  expect_equal(step^2, 5/2)
  expect_true(is.na(interpolationStepK(c(1, 2), c(TRUE, FALSE), 2, .5, 2, 2, noise)))
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
  # Its only available stencil exceeds the local-step cap
  expect_error(step.K(f, 1, cores = 1), "finite derivative")
})

test_that("Kostyrka's method returns reasonable values with tricky functions", {
  s <- step.K(function(x) x^3+1/x, 1)
  # The third derivative vanishes, so a fallback is legitimate here too
  expect_equal(s$value, 2, tolerance = 1e-10)
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

test_that("Kostyrka returns an evaluated step when interpolation is unavailable", {
  f <- function(z) ifelse(abs(z-1) < 1e-10 | abs(z-1) >= 0.15, z^3, NA_real_)
  s <- suppressWarnings(step.K(f, 1, range = c(2^-46, 2^14), cores = 1))
  expect_true(s$par %in% s$iterations$h)
  expect_equal(s$value, (f(1+s$par) - f(1-s$par))/(2*s$par), tolerance = 1e-12)
})

test_that("Kostyrka's fallback is unchanged by different function units", {
  base <- suppressWarnings(step.K(function(z) z^2, 1e-8, range = c(2^-60, 2^14)))
  scaled <- suppressWarnings(step.K(function(z) -2^40*z^2, 1e-8, range = c(2^-60, 2^14)))
  expect_equal(c(scaled$par/base$par, scaled$value/(-2^40*base$value)), c(1, 1))
  expect_null(scaled$iterations$fit)
})

test_that("Kostyrka's fallback fails explicitly without a finite derivative", {
  f <- function(z) ifelse(abs(z) < 1e-160, sign(z)*1e150, NA_real_)
  # Function values are finite on part of the grid, but every available derivative overflows
  expect_error(suppressWarnings(step.K(f, 0, range = c(1e-170, 1), cores = 1)), "finite derivative")
})

test_that("Kostyrka's grid starts at the requested lower step", {
  s <- step.K(sin, 1, h0 = 0.003, shrink.factor = 1/sqrt(2), cores = 1)
  h <- s$iterations$h
  expect_equal(h[1], 0.003 * 2^-36)
  expect_true(max(h) <= s$args$range[2])
})

test_that("Missing estimates cannot join separate truncation ranges", {
  f <- function(z) ifelse(log2(abs(z)) %in% c(-12, -11, -10, -8, -7, -6, -4, -3, -2), z^3, NA_real_)
  s <- suppressWarnings(step.K(f, 0, cores = 1))
  error <- s$iterations$est.error[, "trunc"]
  expect_equal(c(sum(is.finite(error) & error > 0), length(s$args$i.good)), c(6, 0))
})

test_that("An unsuccessful kink fit uses the fallback", {
  skip_if_not_installed("testthat", "3.1.7")
  local_mocked_bindings(optim = function(par, fn, ...) list(par = par, value = fn(par), convergence = 1L), .package = "stats")
  s <- step.K(sin, 1, cores = 1)
  expect_true(!is.null(s$iterations$fallback) && !s$iterations$fit$accepted)
})

test_that("A kink outside the observations uses the fallback", {
  skip_if_not_installed("testthat", "3.1.7")
  local_mocked_bindings(optim = function(par, fn, ...) list(par = c(100, 0), value = 0, convergence = 0L), .package = "stats")
  s <- step.K(sin, 1, cores = 1)
  expect_identical(s$exitcode, 3L)
})


test_that("Parallelisation in Kostyrka's algorithm works", {
  expect_identical(step.K(sin, 1, cores = 1), step.K(sin, 1, cores = 2))
  clus <- parallel::makePSOCKcluster(2)
  expect_identical(step.K(sin, 1, cores = 1), step.K(sin, 1, cl = clus))
  parallel::stopCluster(clus)

  # Testing a slow function
  # f <- function(x) {Sys.sleep(0.1); sin(x)}
  # system.time(step.K(f, 1))
  # system.time(step.K(f, 1, cores = 4))
})
