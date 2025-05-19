test_that("checkDimensions works with good inputs", {
  f <- sin
  good <- c(elementwise = TRUE, vectorised = TRUE, multivalued = FALSE)
  noAttr <- function(x) { # Removing everything except for names
    attributes(x) <- attributes(x)["names"]
    x
  }
  expect_equal(noAttr(checkDimensions(f, 1:3)), good)
  expect_equal(noAttr(checkDimensions(f, 1:3, cl = NULL, cores = 2)), good)

  expect_error(checkDimensions(f, 1:3, elementwise = "yes"), "elementwise")
  expect_error(checkDimensions(f, 1:3, vectorised = "nah"), "vectorised")
  expect_error(checkDimensions(f, 1:3, multivalued = "wat"), "multivalued")

  expect_error(checkDimensions(f, 1:3, multivalued = TRUE, elementwise = TRUE),
               "cannot be elementwise")
  expect_error(checkDimensions(f, 1:3, multivalued = FALSE, elementwise = FALSE,
                               vectorised = TRUE), "cannot be vectorised")

  f <- function(x) stop("Someone left this debug line here")
  expect_error(checkDimensions(f, 1), "cannot be computed")
  expect_warning(checkDimensions(f, 1:3), "at least one finite numeric")

  expect_equal(unname(noAttr(checkDimensions(sin, 1))), c(TRUE, TRUE, FALSE))
  f <- function(x) c(sin(x), cos(x))
  expect_equal(unname(noAttr(checkDimensions(f, 1))), c(FALSE, TRUE, TRUE))
  expect_equal(unname(noAttr(checkDimensions(f, 1:2))), c(FALSE, TRUE, TRUE))
  expect_equal(unname(noAttr(checkDimensions(f, 1:3))), c(FALSE, TRUE, TRUE))

})

test_that("Checks can be skipped without consequences (user responsibility)", {
  expect_equal(unname(checkDimensions(sin, 1:2, elementwise = TRUE, vectorised = TRUE, multivalued = FALSE)),
               c(TRUE, TRUE, FALSE))
  expect_equal(unname(checkDimensions(sin, 1:2, elementwise = FALSE, vectorised = FALSE, multivalued = TRUE)),
               c(FALSE, FALSE, TRUE))
})
