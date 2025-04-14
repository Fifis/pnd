test_that("Jacobians are correct", {
  f <- function(x) c(sine = sum(sin(x)), expon = sum(exp(x)))
  expect_equal(as.numeric(Jacobian(x = 1:3, f)), as.numeric(rbind(sine = cos(1:3), expon = exp(1:3))),
               tolerance = 1e-9)

  x <- structure(1:3, names = c("A", "B", "C"))
  g <- Jacobian(f, x)
  expect_equal(colnames(g), names(x))
  expect_equal(rownames(g), c("sine", "expon"))
  expect_equal(attr(g, "step.size.method"), "default")
  expect_equal(attr(Jacobian(f, x, h = 0.01), "step.size.method"), "user-supplied")
})

test_that("function dimension check works", {
  expect_error(Jacobian(sum, 1:3), "for scalar-valued functions")
})

test_that("Jacobian works with automatic step sizes", {
  f <- function(x) c(x^2 - 2*x + 2, exp(x))
  expect_error(Jacobian(f, x = 0.75, h = "CR"), "step selection works only when")
})

test_that("Jacobian can accept dot arguments", {
  # ... must be accepted and used by checkDimensions() and gradstep()
  f <- function(x, a0) c(sin(x + a0), cos(x + a0))
  f1 <- function(x) c(sin(x + 1), cos(x + 1))
  expect_equal(Jacobian(f, x = 1, a0 = 1, h = 1e-5), Jacobian(f1, x = 1, h = 1e-5))
})


# TODO: compatibility with numDeriv

# TODO: parallelisation
