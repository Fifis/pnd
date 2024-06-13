test_that("Jacobians are correct", {
  f <- function(x) c(sine = sum(sin(x)), expon = sum(exp(x)))
  expect_equal(Jacobian(x = 1:3, f, report = 0), rbind(sine = cos(1:3), expon = exp(1:3)),
               tolerance = 1e-9)

  x <- structure(1:3, names = c("A", "B", "C"))
  g <- Jacobian(f, x)
  expect_equal(colnames(g), names(x))
  expect_equal(rownames(g), c("sine", "expon"))
  expect_equal(attr(g, "step.size.method"), "default")
  expect_equal(attr(Jacobian(f, x, h = 0.01), "step.size.method"), "user-supplied")
})

test_that("function dimension check works", {
  expect_warning(Jacobian(sum, 1:3), "scalar-valued function")
})


# TODO: compatibility with numDeriv

# TODO: parallelisation

