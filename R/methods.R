#' Print method for numerical derivatives
#' @rdname GenD
#' @order 2
#' @export
print.GenD <- function(x) {
  nc <- NCOL(x)
  nr <- NROW(x)
  a <- attributes(x)
  # Printing gradients because those are manageable
  same.h <- all(a$step.size[1] == a$step.size)
  if (nc == 1) {
    cat("Estimated ", if (nr == 1) "derivative" else "gradient", ": [", pasteAnd(printE(x, 4)),
        "].\n(", a$step.size.method, " step size: ",
        if (same.h) printE(a$step.size[1], 1) else pasteAnd(printE(a$step.size, 1)), ").\n", sep = "")
  } else { # Jacobian printed by row
    cat("Estimated Jacobian:\n")
    for (i in 1:nr) cat(paste(printE(x[i, ], 3), collapse = " | "), "\n")
    cat("(", a$step.size.method, " step size: ", if (same.h) printE(a$step.size, 1) else
      paste0(printE(min(a$step.size), 1), "--", printE(max(a$step.size), 1)), ".)\n", sep = "")
  }
}

#' Print method for Hessians
#' @rdname hessian
#' @order 2
#' @export
print.hessian <- function(x) {
  nc <- NCOL(x)
  nr <- NROW(x)
  a <- attributes(x)
  # Printing gradients because those are manageable
  same.h <- all(a$step.size[1] == a$step.size)
  cat("Estimated Hessian:\n")
  for (i in 1:nr) cat(paste(printE(x[i, ], 3), collapse = " | "), "\n")
  cat("Step size: [", if (same.h) printE(a$step.size, 1) else pasteAnd(printE(a$step.size, 1)), "].\n", sep = "")
}

#' Print method for step size
#' @rdname gradstep
#' @order 2
#' @export
#'
#' @examples
#' print(step.CR(x = 1, sin))
#' print(step.DV(x = 1, sin))
#' print(step.plugin(x = 1, sin))
#' print(step.SW(x = 1, sin))
#' print(step.M(x = 1, sin))
#' print(step.K(x = 1, sin))
print.stepsize <- function(x) {
  cat("Step size: ", x$par, " (numerical derivative value: ", x$value, ").\n", sep = "")
  if (x$method == "plug-in") {
    cat(x$counts, " plug-in calculations terminated with code ", x$exitcode, ".\n", sep = "")
  } else if (x$method %in% c("Mathur", "Kink")) {
    cat(x$method, " grid search across ", x$counts, " step sizes ended with code ", x$exitcode, ".\n", sep = "")
  } else {
    cat(x$method, " search terminated after ", paste0(x$counts, collapse = "+"), " iterations with code ", x$exitcode, ".\n", sep = "")
  }
  cat("Error estimates: truncation ", sprintf("%1.0e", x$abs.error[1]),
      ", rounding ", sprintf("%1.0e", x$abs.error[2]), ", total ", sprintf("%1.0e", sum(x$abs.error)), ".\n", sep = "")
}


#' Print method for gradient step size
#' @rdname gradstep
#' @order 2
#' @export
#'
#' @examples
#' f <- function(x) x[1]^3 + sin(x[2])*exp(x[3])
#' print(gradstep(x = c(2, pi/4, 0.5), f))
print.gradstep <- function(x) {
  if (length(x$par) < 2) {
    print.stepsize(x)
  } else {
    all.equal.ss <- all(diff(x$par) == 0)
    cat("Gradient step size: [", if (all.equal.ss) printE(x$par[1], 1) else pasteAnd(printE(x$par, 1)), "].\n",
        "Numerical derivative value: [", pasteAnd(printE(x$value)), "].\n", sep = "")
    all.equal.ec <- all(diff(x$exitcode) == 0)
    if (x$method == "plug-in") {
      cat("Plug-in calculations terminated with code ", if (all.equal.ec) x$exitcode[1] else pasteAnd(x$exitcode), ".\n", sep = "")
    } else if (x$method %in% c("Mathur", "Kink")) {
      cat(x$method, " grid searches across ", x$counts, " step sizes ended with codes ", pasteAnd(x$exitcode), ".\n", sep = "")
    } else {
      cat(x$method, " searches terminated with codes ", pasteAnd(x$exitcode), ".\n", sep = "")
    }
    cat("Estimated errors: truncation [", pasteAnd(sprintf("%1.0e", x$abs.error[, 1])),
        "], rounding [", pasteAnd(sprintf("%1.0e", x$abs.error[, 2])), "].\n", sep = "")
  }
}

#' Print method for check results
#' @rdname checkDimensions
#' @order 2
#' @export
print.checkDimensions <- function(x) {
  cat("Function properties: ", if (!x["elementwise"]) "NOT ", "element-wise, ",
      if (!x["vectorised"]) "NOT ", "vectorised, ",
      if (!x["multivalued"]) "single-valued." else "multi-valued.", sep = "")
}

