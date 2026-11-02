test_that("Jacobians are correct", {
  f <- function(x) c(sine = sum(sin(x)), expon = sum(exp(x)))
  expect_equal(as.numeric(Jacobian(x = 1:3, f)), as.numeric(rbind(sine = cos(1:3), expon = exp(1:3))),
               tolerance = 1e-9)

  x <- structure(1:3, names = c("A", "B", "C"))
  g <- Jacobian(f, x)
  expect_identical(colnames(g), names(x))
  expect_identical(rownames(g), c("sine", "expon"))
  expect_identical(attr(g, "step.size.method"), "default")
  expect_identical(attr(Jacobian(f, x, h = 0.01), "step.size.method"), "user-supplied")
})

test_that("Jacobian works for scalar functions", {
  h1 <- function(x) {
    h <- rep(NA, 1)
    h[1] <- x[1] + x[2] + x[3] - 1
    h
  }
  h2 <- function(x) {
    h <- rep(NA, 1)
    h[1] <- 6*x[2] + 4*x[3] - x[1]^3 - 3
    h[2] <- x[1]
    h[3] <- x[2]
    h[4] <- x[3]
    h
  }
  p0 <- structure(c(0.1, 0.2, 0.3), names = c("A", "B", "C"))
  h1.jac <- Jacobian(h1, p0)
  h2.jac <- Jacobian(h2, p0)
  expect_true(inherits(h1.jac, "GenD"))
  expect_true(inherits(h2.jac, "GenD"))
  expect_identical(dim(h1.jac), c(1L, 3L))
  expect_identical(dim(h2.jac), c(4L, 3L))
  expect_identical(colnames(h1.jac), c("A", "B", "C"))
  expect_identical(colnames(h2.jac), c("A", "B", "C"))
})

test_that("Jacobian works well with multi-valued functions", {
  set.seed(1)
  X <- rchisq(20, 3)
  Y <- 1 + X + rnorm(20)
  d <- data.frame(Y, X)
  mod <- lm(Y ~ X, data = d)
  # Least-squares moment conditions
  g <- function(theta, data) {
    res <- Y - theta[1] - theta[2]*X
    cbind(U = res, UX = res*data$X)
  }
  # colMeans(g(mod$coefficients, d))  # Around 1e-16
  ghat <- function(theta) colMeans(g(theta, data = d))
  # ghat(mod$coefficients)
  jac <- Jacobian(ghat, x = mod$coefficients, elementwise = FALSE, multivalued = TRUE, vectorised = FALSE)
  expect_equal(-c(1, mean(X), mean(X), mean(X^2)), as.numeric(jac), tolerance = 1e-8)
  # TODO: make it work without checks
})

test_that("Jacobian works with automatic step sizes", {
  f <- function(x) c(x^2 - 2*x + 2, exp(x))
  expect_error(Jacobian(f, x = 0.75, h = "CR"), "step selection works only when")
})

test_that("Jacobian can accept dot arguments", {
  # ... must be accepted and used by checkDimensions() and gradstep()
  f <- function(x, a0) c(sin(x + a0), cos(x + a0))
  f1 <- function(x) c(sin(x + 1), cos(x + 1))
  expect_identical(Jacobian(f, x = 1, a0 = 1, h = 1e-5), Jacobian(f1, x = 1, h = 1e-5))
})

test_that("Jacobian can work on an arbitrary stencil", {
  f <- function(x) c(sum(sin(x[1:2])), exp(x[3]))
  jac <- rbind(c(cos(1:2), 0), c(0, 0, exp(3)))
  d2 <- Jacobian(f, 1:3, h = 1e-4, stencil = c(-1, 1))
  d4 <- Jacobian(f, 1:3, h = 1e-4, stencil = c(-2, -1, 1, 2))  # More accurate for a large step
  expect_equal(d2, d4, tolerance = 1e-7)
  # Comparing errors of non-zero elements
  dd2 <- abs(d2-jac)
  dd4 <- abs(d4-jac)
  zeros <- jac == 0
  expect_true(all(dd2[!zeros] > dd4[!zeros]))
})

test_that("Jacobians correctly handle ambiguous cases Rn -> Rn", {
  f <- function(x) c(sin(x[1]) + cos(x[2]), exp(x[1]*x[2]))
  x0 <- c(1, 1.5)
  d2 <- matrix(as.numeric(Jacobian(f, x0)), nrow = 2)
  jac <- matrix(c(cos(1), 1.5*exp(1.5), -sin(1.5), exp(1.5)), nrow = 2)
  expect_equal(d2, jac, tolerance = 1e-8)
})

# Tombs of the ancients: numDeriv compatibility

test_that("Jacobian resolves func before checking square mappings", {
  f <- function(x, scale0) c(total = scale0*sum(x), product = scale0*prod(x))
  x <- 1:2
  expect_warning(d <- Jacobian(func = f, x = x, scale0 = 3, h = 1/32), "Use the argument 'FUN'")
  expect_equal(as.numeric(d), c(3, 6, 3, 3), tolerance = 1e-10)
  expect_identical(d, Jacobian(FUN = f, x = x, scale0 = 3, h = 1/32))
})

test_that("Jacobian preserves names and supplied values with func", {
  f <- function(x, scale0) c(total = scale0*sum(x), product = scale0*prod(x))
  x <- c(A = 1, B = 2)
  direct <- Jacobian(FUN = f, x = x, scale0 = 3, h = 1/32)
  for (f0 in list(NULL, f(x, scale0 = 3))) {
    expect_warning(d <- Jacobian(func = f, x = x, scale0 = 3, h = 1/32, f0 = f0), "Use the argument 'FUN'")
    expect_identical(d, direct)
  }
})

test_that("Jacobian validates function arguments and gives FUN precedence", {
  f <- function(x) c(sum(x), prod(x))
  unused <- function(x) stop("The compatibility alias should not be used.")
  expect_identical(Jacobian(FUN = f, func = unused, x = 1:2, h = 1/32), Jacobian(FUN = f, x = 1:2, h = 1/32))
  expect_error(Jacobian(x = 1:2), "Pass the function")
  expect_error(Jacobian(func = 1, x = 1:2), "Pass the function")
  expect_error(Jacobian(FUN = 1, func = f, x = 1:2), "must be a function")
})

# TODO: parallelisation
