test_that("GenD prints named and unnamed gradients as it should", {
  x <- structure(1:3, names = c("Andy", "Bradley", "Ca"))

  o <- capture.output(print(Grad(function(x) prod(sin(x)), 1)))
  expect_true(grepl("Estimated derivative", o[1]))
  expect_true(grepl("default step", o[2]))

  o <- capture.output(print(Grad(function(x) prod(sin(x)), x)))
  expect_true(grepl("Estimated gradient", o[1]))
  expect_true(grepl("Bradley", o[2]))

  o <- capture.output(print(Jacobian(function(x) c(prod(sin(x)), sum(exp(x))), x)))
  expect_true(grepl("step size range", o[length(o)]))
})

test_that("Hessian rudimentary print support", {
  o <- capture.output(print(Hessian(function(x) prod(sin(x)), 1:3)))
  expect_true(grepl("Estimated Hessian", o[1]))
  expect_true(grepl("step size", o[length(o)]))
})

test_that("Step-size search prints uniform and adequate results", {
  o1 <- capture.output(print(step.CR(x = 1, sin)))
  o2 <- capture.output(print(step.DV(x = 1, sin)))
  o3 <- capture.output(print(step.plugin(x = 1, sin)))
  o4 <- capture.output(print(step.SW(x = 1, sin)))
  o5 <- capture.output(print(step.M(x = 1, sin)))
  o6 <- capture.output(print(step.K(x = 1, sin)))
  o <- list(o1, o2, o3, o4, o5, o6)
  expect_true(all(grepl("Step size:", sapply(o, "[", 1))))
  expect_true(all(grepl("(search|calculations) (terminated|across)", sapply(o, "[", 2))))
  expect_true(all(grepl("Error estimates:", sapply(o, "[", 3))))

  f <- function(x) x[1]^3 + sin(x[2])*exp(x[3])
  o <- capture.output(print(gradstep(x = c(2, pi/4, 0.5), f)))
  expect_true(grepl("Gradient step size", o[1]))
  expect_true(grepl("terminated with codes", o[3]))
  expect_true(grepl("Estimated errors", o[4]))
})

test_that("Dimension checks are printed adequately", {
  d <- capture.output(print(checkDimensions(sin, x = 1)))
  expect_equal(d, "Function properties: element-wise, vectorised, single-valued.")
})


