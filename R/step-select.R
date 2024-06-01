#' Curtis--Reid automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of .Machine$double.eps^(1/3)
#'   (or absolute step if \code{x} is zero).
#' @param version Character scalar: \code{"original"} for the original 1974 version by
#'   Curtis and Reid; \code{"modified"} for Kostyrka’s 2024 modification, which adds an
#'   extra evaluation for a more accurate estimate of the truncation error.
#' @param aim Positive real scalar: desired ratio of truncation-to-rounding error. The \code{"original"}
#'   version over-estimates the truncation error, hence a higher \code{aim} is recommended.
#'   For the \code{"modified"} version, aim should be close to 1.
#' @param acc.order Numeric scalar: in the modified version, allows returning a 4th-order accurate
#'   central difference using a 4-point stencil from the last iteration with
#'   half the optimal step size.
#' @param tol Numeric scalar greater than 1: tolerance multiplier for determining when to stop
#'   the algorithm based on the current estimate being between `aim/tol` and `aim*tol`.
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param maxit Integer: maximum number of algorithm iterations to prevent infinite
#'   loops in degenerate cases.
#' @param seq.tol Numeric scalar: maximum relative difference between old and new
#'   step sizes for declaring convergence.
#' @param diagnostics Logical: if \code{TRUE}, returns the full iteration history
#'   including all function evaluations.
#' @param ... Passed to \code{FUN}.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{curtis1974choice}{pnd} algorithm.
#' If the estimated third derivative is exactly zero, then, the initial step size
#' is multipled by 2 and returned.
#'
#' TODO: mention that f must be one-dimensional
#'
#' @return A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
#'   step size found, \code{value} -- the estimated numerical first derivative (central differences;
#'   very useful for computationally expensive functions), \code{counts} -- the number of
#'   iterations (each iteration includes three function evaluations), \code{exitcode} --
#'   an integer code indicating the termination status:
#'   \code{0} indicates optimal termination within tolerance,
#'   \code{1} means that the third derivative is zero (large step size preferred),
#'   \code{2} is returned if there is no change in step size within tolerance,
#'   \code{3} indicates a solution at the boundary of the allowed value range,
#'   \code{4} signals that the maximum number of iterations was reached.
#'   \code{message} -- a summary message of the exit status.
#'   If \code{diagnostics} is \code{TRUE}, \code{iterations} is a list
#'   including the full step size search path, argument grids, function values on
#'   those grids, estimated error ratios, and estimated derivative values.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.CR(x = 2, f)
#' step.CR(x = 2, f, h0 = 1e-3, diagnostics = TRUE)
#' step.CR(x = 2, f, version = "modified")
step.CR <- function(x, FUN, h0 = 6e-6 * (abs(x) + (x == 0)),
                    version = c("original", "modified"),
                    aim = if (version[1] == "original") 100 else 1,
                    acc.order = c(2L, 4L),
                    tol = if (version[1] == "original") 10 else 4,
                    range = NULL, maxit = 20L, seq.tol = 1e-8,
                    diagnostics = FALSE, ...) {
  version <- version[1]
  if (!(version %in% c("original", "modified"))) stop("step.CR: 'version' must be either 'original' or 'modified'.")
  vanilla <- version == "original"
  if (is.null(range)) range <- h0 / c(1e5, 1e-5)
  if (range[1] < 2*.Machine$double.eps) range[1] <- 2*.Machine$double.eps
  if (tol <= 1) stop("The tolerance must be a positive number greater than 1 (e.g. 4).")
  acc.order <- acc.order[1]
  if (acc.order == 4 && vanilla) {
    warning("The original 1974 algorithm does not support 4th-order accuracy. Setting acc.order = 2.")
    acc.order <- 2
  }
  if (acc.order == 4) aim <- aim * .Machine$double.eps^(-2/15)
  target <- sort(c(aim / tol, aim * tol))

  getRatio <- function(FUN, x, h, vanilla, ...) {
    # TODO: parallelise
    xgrid <- x + if (vanilla) c(-h, 0, h) else c(-h, -h/2, h/2, h)
    fgrid  <- sapply(xgrid, FUN, ...)
    if (vanilla) {
      fd <- (fgrid[3] - fgrid[2]) / h
      bd <- (fgrid[2] - fgrid[1]) / h
      cd <- (fgrid[3] - fgrid[1]) / h * 0.5
      etrunc <- abs(cd - fd)
    } else {
      cd <- (fgrid[4] - fgrid[1]) / h * 0.5
      cd.half <- (fgrid[3] - fgrid[2]) / h
      cd4 <- sum(fgrid * c(1, -8, 8, -1)) / h / 6
      etrunc <- abs(cd - cd.half) * 4 / 3
    }
    eround <- 0.5 * max(abs(fgrid)) * .Machine$double.eps / h
    u <- etrunc / eround
    deriv <- if (vanilla) c(cd = cd, bd = bd, fd = fd) else c(cd = cd, cd.half = cd.half, cd4 = cd4)
    ret <- list(h = h, x = xgrid, f = fgrid, ratio = u, deriv = deriv,
                est.error = c(trunc = etrunc, round = eround))
    return(ret)
  }

  ulist <- list()
  i <- 1
  exitcode <- 0

  while (i <= maxit) {
    hold <- if (i > 1) hnew else NA
    hnew <- if (i > 1) hold * (aim / max(ulist[[i-1]]$ratio, 1))^(1/3) else h0

    # Check the relative change of the step size, which is possible at
    # the 2nd iteration even before the 2nd error calculation
    if (i > 1) {
      if (max(hnew/hold, hold/hnew) - 1 < seq.tol) {
        exitcode <- 2
        break # Step 4: if the step size does not change, stop
      } else { # 4a, b: outside the range, replace with the border
        if (hnew <= range[1]) hnew <- range[1]
        if (hnew >= range[2]) hnew <- range[2]
        if (max(hnew/hold, hold/hnew) - 1 < seq.tol) {
          exitcode <- 3
          break
        }
      }
    }

    ulist[[i]] <- getRatio(FUN = FUN, x = x, h = hnew, vanilla = vanilla, ...)
    if (any(bad <- !is.finite(ulist[[i]]$fgrid))) {
      stop(paste0("Could not compute the function value at ", paste0(ulist[[i]]$xgrid[bad], collapse = " and "),
                  ". Change the range, which is currently [", paste0(range, collapse = "; "),
                  "], and/or try a different starting h0, which is currently ", h0, "."))
    }

    if (ulist[[i]]$ratio >= target[1] && ulist[[i]]$ratio <= target[2]) {
      break # Successful termination by matching the range
    }
    if (ulist[[i]]$ratio == 0) { # Zero truncation error -> only rounding error
      exitcode <- 1
      hnew <- hnew * 4 # For linear or quadratic functions, a large h is preferred
      break
    }

    i <- i+1
  }
  i <- length(ulist)
  if (i >= maxit) {
    exitcode <- 4
  }
  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation error is exactly zero, large step is favoured",
                "step size did not change between iterations",
                paste0("step size landed on the range ", if (hnew == range[1]) "left" else
                       "right", " end; consider extending the range"),
                "maximum number of iterations reached")
  if (diagnostics) {
    diag.list <- list(h = do.call(c, lapply(ulist, "[[", "h")),
                      x = do.call(rbind, lapply(ulist, "[[", "x")),
                      f = do.call(rbind, lapply(ulist, "[[", "x")),
                      ratio = do.call(c, lapply(ulist, "[[", "ratio")),
                      deriv = do.call(rbind, lapply(ulist, "[[", "deriv")))
  }
  ret <- list(par = hnew,
              value = if (acc.order == 4) unname(ulist[[i]]$deriv["cd4"]) else unname(ulist[[i]]$deriv["cd"]),
              counts = i, exitcode = exitcode, message = msg,
              iterations = if (diagnostics) diag.list else NULL)
  return(ret)
}


#' Dumontet--Vignes automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of .Machine$double.eps^(1/3)
#'   (or absolute step if \code{x} is zero).
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param alpha Numeric scalar >= 1 indicating the relative reduction in the
#'   number of accurate bits due to the calculation of \code{FUN}. A value of \code{1}
#'   implies that \code{FUN(x)} is assumed to have all bits accurate with maximum relative
#'   error of \code{.Machine$double.eps/2}. A value of \code{2} indicates that the number of
#'   accurate bits is half the mantissa length, \code{4} if it is quarter etc. The algorithm authors
#'   recommend \code{4/3} even for highly accurate functions.
#' @param ratio.limits Numeric vector of length 4 defining the acceptable ranges
#'   for step size: the algorithm stops if the relative perturbation of the third derivative by
#'   amplified rounding errors falls either between the 1st and 2nd elements or between
#'   the 3rd and 4th elements.
#' @param maxit Maximum number of algorithm iterations to avoid infinite loops in cases
#'   the desired relative perturbation factor cannot be achieved within the given \code{range}.
#'   Consider extending the range if this limit is reached.
#' @param diagnostics Logical: if \code{TRUE}, returns the full iteration history
#'   including all function evaluations.
#'   Note: the history tracks the third derivative, not the first.
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{dumontet1977determination}{pnd} algorithm.
#' If the estimated third derivative is exactly zero, the function assumes a third
#' derivative of 1 to prevent division-by-zero errors.
#'
#'
#' @return A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
#'   step size found, \code{value} -- the estimated numerical first derivative (central
#'   differences), \code{counts} -- the number of iterations (each iteration includes
#'   four function evaluations), \code{exitcode} -- an integer code indicating the
#'   termination status:
#'   \code{0} indicates optimal termination within tolerance,
#'   \code{1} means that the third derivative is zero (large step size preferred),
#'   \code{3} indicates a solution at the boundary of the allowed value range,
#'   \code{4} signals that the maximum number of iterations was reached.
#'   \code{message} is a summary message of the exit status.
#'   If \code{diagnostics} is \code{TRUE}, \code{iterations} is a list
#'   including the full step size search path (NB: for the 3rd derivative),
#'   argument grids, function values on those grids, and estimated 3rd derivative values.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.DV(x = 2, f)
#' step.DV(x = 2, f, h0 = 1e-3, diagnostics = TRUE)
step.DV <- function(x, FUN, h0 = 6e-6 * (abs(x) + (x == 0)),
                    range = NULL, alpha = 4/3,
                    ratio.limits = c(1/15, 1/2, 2, 15),
                    maxit = 40L, diagnostics = FALSE, ...) {
  k0 <- h0 * .Machine$double.eps^(-2/15)
  if (is.null(range)) range <- h0 / c(1e6, 1e-6)
  range3 <- range * .Machine$double.eps^(-2/15)
  P <- 2^(log2(.Machine$double.eps/2) / alpha)

  getRatio <- function(FUN, x, k, P, ...) {
    # TODO: parallelise
    xgrid <- x + c(-2*k, -k, k, 2*k)
    fgrid  <- sapply(xgrid, FUN, ...)
    tgrid <- fgrid * c(-0.5, 1, -1, 0.5) # T1, ..., T4 from the paper
    A <- sum(tgrid[tgrid > 0]) # Only positive terms
    B <- sum(tgrid[tgrid < 0]) # Only negative terms
    # TODO: error if fgrid has different sign
    f3inf <- (A / (1+P) + B / (1-P)) / k^3
    f3sup <- (A / (1-P) + B / (1+P)) / k^3
    f3 <- (f3inf + f3sup) / 2 # Estimate of third derivative
    L <- f3sup / f3inf
    ret <- list(k = k, x = xgrid, f = fgrid, ratio = L,
                deriv = c(f3inf = f3inf, f3 = f3, f3sup = f3sup))
    return(ret)
  }

  ulist <- list()
  i <- 1
  exitcode <- 0
  ndownwards <- 0 # For counting the number of downwards shrinkages

  while (i <= maxit) {
    if (i == 1) k <- k0
    ulist[[i]] <- getRatio(FUN = FUN, x = x, k = k, P = P, ...)
    if (any(bad <- !is.finite(ulist[[i]]$fgrid))) {
      stop(paste0("Could not compute the function value at ", paste0(ulist[[i]]$xgrid[bad], collapse = " and "),
                  ". Change the range, which is currently [", paste0(range, collapse = "; "),
                  "], and/or try a different starting h0, which is currently ", h0, "."))
    }

    if (ulist[[i]]$ratio < ratio.limits[1] || ulist[[i]]$ratio > ratio.limits[4]) {
      # The numerical error is too high
      range3[1] <- k
    } else {
      if (ulist[[i]]$ratio > ratio.limits[2] && ulist[[i]]$ratio < ratio.limits[3]) {
        # The numerical error is too small
        range3[2] <- k
        ndownwards <- ndownwards + 1
      } else {
        # The numerical error is in the good range, the step size is optimal
        break
      }
    }
    k <- sqrt(prod(range3)) # The interval has been sub-divided; take its geometric mean
    i <- i+1
  }

  i <- length(ulist)
  if (i >= maxit) {
    # Was the procedure systematically unsuccsessul?
    exitcode <- if (maxit > 5 && (ndownwards >= maxit-2 || ndownwards <= 2)) 3 else 4
  }

  f3 <- sum(ulist[[i]]$f * c(-0.5, 1, -1, 0.5)) / k^3
  if (f3 == 0) {
    exitcode <- 1
    f3 <- 1
  }

  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation error is exactly zero, large step is favoured",
                "",
                paste0("step size too close to the range ", if (ndownwards > maxit/2) "left" else
                       "right", " end; consider extending the range"),
                "maximum number of iterations reached")

  f0 <- mean(ulist[[i]]$f[2:3]) # Approximately f(x); the error is minuscule and immaterial
  h <- (1.67 * P * abs(f0/f3))^(1/3) # Formula 36 from Dumontet & Vignes (1977)
  h <- h + x # Minimising the representation error
  h <- h - x
  xgrid <- c(x-h, x+h)
  fgrid <- sapply(xgrid, FUN, ...)
  cd <- (fgrid[2] - fgrid[1]) / h / 2
  err <- P*abs(f0)/h/3 + (exitcode != 1) * (h^5*f3^2/36/P/abs(f0) - h^8*abs(f3)^3/648/P^2/f0^2)

  if (diagnostics) {
    diag.list <- list(k = do.call(c, lapply(ulist, "[[", "k")),
                      x = do.call(rbind, lapply(ulist, "[[", "x")),
                      f = do.call(rbind, lapply(ulist, "[[", "x")),
                      ratio = do.call(c, lapply(ulist, "[[", "ratio")),
                      deriv3 = do.call(rbind, lapply(ulist, "[[", "deriv")))
  }

  ret <- list(par = h, value = cd,
              counts = i, exitcode = exitcode, message = msg,
              iterations = if (diagnostics) diag.list else NULL,
              abs.error = err)
  return(ret)

}
