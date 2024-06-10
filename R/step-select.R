.safeF <- function(FUN, x, ...) tryCatch(FUN(x, ...), error = function(e) return(NA))
.pasteAnd <- function(x) paste0(x, collapse = ", ")

#' Curtis--Reid automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param version Character scalar: \code{"original"} for the original 1974 version by
#'   Curtis and Reid; \code{"modified"} for Kostyrka’s 2024 modification, which adds an
#'   extra evaluation for a more accurate estimate of the truncation error.
#' @param aim Positive real scalar: desired ratio of truncation-to-rounding error. The \code{"original"}
#'   version over-estimates the truncation error, hence a higher \code{aim} is recommended.
#'   For the \code{"modified"} version, aim should be close to 1.
#' @param acc.order Numeric scalar: in the modified version, allows searching for a
#'   step size that would be optimal for a 4th-order-accurate central difference
#'   See the Details section below.
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
#' is multiplied by 4 and returned.
#'
#' If 4th-order accuracy (4OA) is requested, then, two things happen. Firstly,
#' since 4OA differences requires a larger step size and the truncation error for
#' the 2OA differences grows if the step size is larger than the optimal one,
#' a higher ratio of truncation-to-rounding errors should be targeted. Secondly,
#' a 4OA numerical derivative is returned, but the truncation and rounding errors
#' are still estimated for the 2OA differences. Therefore, the estimating truncation
#' error is higher and the real truncation error of 4OA differences is lower.
#'
#' TODO: mention that f must be one-dimensional
#'
#' @return A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
#'   step size found, \code{value} -- the estimated numerical first derivative (central differences;
#'   very useful for computationally expensive functions), \code{counts} -- the number of
#'   iterations (each iteration includes three function evaluations), \code{abs.error} --
#'   an estimate of the total approximation error (sum of truncation and rounding errors),
#'   \code{exitcode} -- an integer code indicating the termination status:
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
#' step.CR(x = 2, f, version = "modified", acc.order = 4)
step.CR <- function(x, FUN, h0 = 1e-5 * (abs(x) + (x == 0)),
                    version = c("original", "modified"),
                    aim = if (version[1] == "original") 100 else 1,
                    acc.order = c(2L, 4L),
                    tol = if (version[1] == "original") 10 else 4,
                    range = h0 / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                    diagnostics = FALSE, ...) {
  version <- version[1]
  if (!(version %in% c("original", "modified"))) stop("step.CR: 'version' must be either 'original' or 'modified'.")
  vanilla <- version == "original"
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
    fgrid  <- sapply(xgrid, function(z) .safeF(FUN, z, ...))
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
      if (abs(hnew/hold - 1) < seq.tol) {
        exitcode <- 2
        break # Step 4: if the step size does not change, stop
      } else { # 4a, b: outside the range, replace with the border
        if (hnew < range[1]) hnew <- range[1]
        if (hnew > range[2]) hnew <- range[2]
        if (max(hnew/hold, hold/hnew) - 1 < seq.tol) {
          exitcode <- 3
          break
        }
      }
    }

    ulist[[i]] <- getRatio(FUN = FUN, x = x, h = hnew, vanilla = vanilla, ...)
    if (any(bad <- !is.finite(ulist[[i]]$fgrid))) {
      stop(paste0("Could not compute the function value at ", .pasteAnd(ulist[[i]]$xgrid[bad]),
                  ". Change the range, which is currently [", .pasteAnd(range),
                  "], and/or try a different starting h0, which is currently ", h0, "."))
    }

    if (ulist[[i]]$ratio >= target[1] && ulist[[i]]$ratio <= target[2]) {
      break # Successful termination by matching the range
    }
    if (ulist[[i]]$ratio == 0) { # Zero truncation error -> only rounding error
      exitcode <- 1
      hnew <- hnew * 4 # For linear or quadratic functions, a large h is preferred
      ulist[[i+1]] <- getRatio(FUN = FUN, x = x, h = hnew, vanilla = vanilla, ...)
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
                paste0("step size landed on the range ", if (ulist[[i]]$h == range[1]) "left" else
                       "right", " end; consider extending the range"),
                "maximum number of iterations reached")
  if (diagnostics) {
    diag.list <- list(h = do.call(c, lapply(ulist, "[[", "h")),
                      x = do.call(rbind, lapply(ulist, "[[", "x")),
                      f = do.call(rbind, lapply(ulist, "[[", "f")),
                      deriv = do.call(rbind, lapply(ulist, "[[", "deriv")),
                      est.error = do.call(rbind, lapply(ulist, "[[", "est.error")),
                      ratio = do.call(c, lapply(ulist, "[[", "ratio")))
  }
  ret <- list(par = ulist[[i]]$h,
              value = if (acc.order == 4) unname(ulist[[i]]$deriv["cd4"]) else unname(ulist[[i]]$deriv["cd"]),
              counts = i, exitcode = exitcode, message = msg,
              abs.error = sum(ulist[[i]]$est.error),
              iterations = if (diagnostics) diag.list else NULL)
  return(ret)
}


#' Dumontet--Vignes automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
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
#'   four function evaluations), \code{abs.error} -- an estimate of the total
#'   approximation error (sum of truncation and rounding errors),
#'   \code{exitcode} -- an integer code indicating the termination status:
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
step.DV <- function(x, FUN, h0 = 1e-5 * (abs(x) + (x == 0)),
                    range = h0 / c(1e6, 1e-6), alpha = 4/3,
                    ratio.limits = c(1/15, 1/2, 2, 15),
                    maxit = 40L, diagnostics = FALSE, ...) {
  k0 <- h0 * .Machine$double.eps^(-2/15)
  range3 <- range * .Machine$double.eps^(-2/15)
  P <- 2^(log2(.Machine$double.eps/2) / alpha)

  getRatio <- function(FUN, x, k, P, ...) {
    # TODO: parallelise
    xgrid <- x + c(-2*k, -k, k, 2*k)
    fgrid  <- sapply(xgrid, function(z) tryCatch(FUN(z, ...), error = function(e) return(NA)))
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
      stop(paste0("Could not compute the function value at ", .pasteAnd(ulist[[i]]$xgrid[bad]),
                  ". Change the range, which is currently [", .pasteAnd(range),
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

  f0 <- mean(ulist[[i]]$f[2:3]) # Approximately f(x); the error is small for small h
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
                      f = do.call(rbind, lapply(ulist, "[[", "f")),
                      ratio = do.call(c, lapply(ulist, "[[", "ratio")),
                      deriv3 = do.call(rbind, lapply(ulist, "[[", "deriv")))
  }

  ret <- list(par = h, value = cd, counts = i, exitcode = exitcode,
              message = msg, abs.error = err,
              iterations = if (diagnostics) diag.list else NULL)
  return(ret)

}


#' Stepleman--Winarsky automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param shrink.factor A scalar greater than 1 that is used to divide the step size
#'   during the search. The authors recommend 4, but this may be result in earlier
#'   termination at slightly sub-optimal steps. Change to 3 for a faster search.
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param seq.tol Numeric scalar: maximum relative difference between old and new
#'   step sizes for declaring convergence.
#' @param max.rel.error Positive numeric scalar > 0 indicating the maximum relative
#'   error of function evaluation. For highly accurate functions with all accurate bits
#'   is equal to half of machine epsilon. For noisy functions (derivatives, integrals,
#'   output of optimisation routines etc.), it is higher.
#' @param maxit Maximum number of algorithm iterations to avoid infinite loops.
#'   Consider trying some smaller or larger initial step size \code{h0}
#'   if this limit is reached.
#' @param diagnostics Logical: if \code{TRUE}, returns the full iteration history
#'   including all function evaluations.
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{stepleman1979adaptive}{pnd} algorithm.
#'
#'
#' @return A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
#'   step size found, \code{value} -- the estimated numerical first derivative (central
#'   differences), \code{counts} -- the number of iterations (each iteration includes
#'   four function evaluations), \code{abs.error} -- an estimate of the total
#'   approximation error (sum of truncation and rounding errors),
#'   \code{exitcode} -- an integer code indicating the termination status:
#'   \code{0} indicates optimal termination within tolerance,
#'   \code{2} is returned if there is no change in step size within tolerance,
#'   \code{3} indicates a solution at the boundary of the allowed value range,
#'   \code{4} signals that the maximum number of iterations was reached.
#'   \code{message} is a summary message of the exit status.
#'   If \code{diagnostics} is \code{TRUE}, \code{iterations} is a list
#'   including the full step size search path,
#'   argument grids, function values on those grids, estimated derivative values,
#'   estimated error values, and monotonicity check results.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4  # The derivative at 1 is 4
#' step.SW(x = 1, f)
#' step.SW(x = 1, f, h0 = 1e-9, diagnostics = TRUE) # Starting too low
#' # Starting somewhat high leads to too many preliminatry iterations
#' step.SW(x = 1, f, h0 = 10, diagnostics = TRUE)
#' step.SW(x = 1, f, h0 = 1000, diagnostics = TRUE) # Starting absurdly high
#'
#' f <- sin  # The derivative at pi/4 is sqrt(2)/2
#' step.SW(x = pi/4, f)
#' step.SW(x = pi/4, f, h0 = 1e-9, diagnostics = TRUE) # Starting too low
#' step.SW(x = pi/4, f, h0 = 0.1, diagnostics = TRUE) # Starting slightly high
#' # The following two example fail because the truncation error estimate is invalid
#' step.SW(x = pi/4, f, h0 = 10, diagnostics = TRUE)   # Warning
#' step.SW(x = pi/4, f, h0 = 1000, diagnostics = TRUE) # Warning
step.SW <- function(x, FUN, h0 = 1e-5 * (abs(x) + (x == 0)),
                    shrink.factor = 2, range = h0 / c(1e12, 1e-8),
                    seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2,
                    maxit = 40L, diagnostics = FALSE, ...) {
  i <- 1

  getRatio <- function(FUN, x, h, do.f0 = FALSE, ratio.last = NULL,
                       ratio.beforelast = NULL, ...) {
    # TODO: parallelise
    xgrid <- x + if (do.f0) c(-h, 0, h) else c(-h, h)
    fgrid  <- sapply(xgrid, function(z) .safeF(FUN, z, ...))
    cd <- (fgrid[length(fgrid)] - fgrid[1]) / h * 0.5
    etrunc <- if (is.null(ratio.last)) NA else # Richardson extrapolation
      abs((cd - ratio.last$deriv) / (1 - (ratio.last$h / h)^2))
    eround <- max.rel.error * max(abs(fgrid)) / h
    dlast <- if (!is.null(ratio.last) && !is.null(ratio.beforelast))
      c(ratio.beforelast$deriv, ratio.last$deriv, cd) else NULL
    monotone.fp  <- if (!is.null(dlast)) prod(diff(dlast)) > 0 else NA
    monotone.dfp <- if (!is.null(dlast))
      abs(dlast[2] - dlast[3]) <= abs(dlast[2] - dlast[1]) else NA
    ret <- list(h = h, x = if (do.f0) xgrid[c(1, 3)] else xgrid,
                f = if (do.f0) fgrid[c(1, 3)] else fgrid,
                deriv = cd, f0 = if (do.f0) fgrid[2] else NULL,
                est.error = c(trunc = etrunc, round = eround),
                monotone = c(monotone.fp, monotone.dfp))
    return(ret)
  }

  exitcode <- 0
  main.loop <- close.left <- FALSE
  first.main <- TRUE
  ulist <- list()

  while (i <= maxit) {
    if (!main.loop) {
      if (i == 1) {
        ulist[[i]] <- getRatio(FUN, x, h0, do.f0 = TRUE, ratio.last = NULL, ...)
        f0 <- ulist[[i]]$f0
        hnew <- ulist[[i]]$h
      }
      if (any(bad <- !is.finite(ulist[[i]]$f))) {
        if (is.na(f0)) {
          stop(paste0("Could not compute the function value at ", x, ". FUN(x) must be finite."))
        } else {
          stop(paste0("Could not compute the function value at ", .pasteAnd(ulist[[i]]$x[bad]),
                      ". FUN(", x, ") is finite -- reduce the step h0, which is currently ", h0, "."))
        }
      }

      # First check: are the function values of different signs?
      # If yes, h is too large -- jump to the step-shrinking main loop:
      if (prod(sign(ulist[[i]]$f)) < 0) { # f(x+h) should be close to f(x-h)
        main.loop <- TRUE
        i.prelim <- i
        next
      } else { # The bisection part
        while (!main.loop && i <= maxit) {
          hold <- hnew
          # N: approx. number of inaccurate digits due to rounding at h0 and hopt
          Nh0 <- log10(abs(f0)) - log10(2) - log10(hnew) - log10(abs(ulist[[i]]$deriv))
          Nhopt <- -log10(max.rel.error) / 3
          rounding.nottoosmall <- Nh0 > 0 # Formula 3.16: some digits were lost,
          # the rounding error not zero, the step size not extremely large
          rounding.small <- Nh0 <= Nhopt - log10(shrink.factor) # Formula 3.15:
          # the rounding error is not worse than the truncation error
          if (rounding.small && rounding.nottoosmall) {
            main.loop <- TRUE # We are in the truncation branch -- go to shrinkage directly
            break
          } else { # Either rounding.small (too small) or rounding.nottoosmall (too large)
            i <- i + 1
            # This is what is necessary for 3.15 to hold with a small safety margin
            hnew <- hnew * 10^(Nh0 - (Nhopt - log10(shrink.factor)) + log10(2))
            # If 3.15 does not hold, this increases the step size
            # If 3.16 does not hold, shrink hnew towards the 3.15-optimal value
            if (!rounding.nottoosmall) hnew <- (hold + hnew) / 2 # Bisection
            if (hnew > range[2]) hnew <- range[2] # Safeguarding against huge steps
            ulist[[i]] <- getRatio(FUN, x, hnew, do.f0 = FALSE, ratio.last = ulist[[i-1]], ...)
            bad <- !is.finite(ulist[[i]]$f)
            if (any(bad) && !rounding.small)
              stop(paste0("Could not compute the function value at ", .pasteAnd(ulist[[i]]$x[bad]),
                          ". FUN(", x, ") is finite -- try a step h0 larger than ",
                          h0, " but smaller than ", hold, " (currently ", hnew, ")."))
            if (any(bad) && !rounding.nottoosmall)
              stop(paste0("Could not compute the function value at ", .pasteAnd(ulist[[i]]$x[bad]),
                          ". FUN(", x, ") is finite -- try a step h0 smaller than ", hnew, "."))
          }
        } # End initial step search
      }
      i.prelim <- i # Passing to the first iteration of the main loop with 2 function values saved
    } # End preliminary loop


    if (first.main) { # First main loop: extra h1 and h2 needed
      for (j in 1:2) { # Try a decreasing sequence
        i <- i + 1
        hnew <- hnew / shrink.factor
        ulist[[i]] <- getRatio(FUN, x, hnew, ratio.last = ulist[[i-1]],
                               ratio.beforelast = if (j == 1) NULL else ulist[[i-2]])
      }
      first.main <- FALSE
    } else { # Monotonicity satisfied, continuing shrinking, only h[i+1] needed
      if (abs(ulist[[i]]$h/ulist[[i-1]]$h - 1) < seq.tol) {
        exitcode <- 2 # If h did not shrink, it must have hit the lower bound
        break  # or something else went wrong; this code will most likely
        # be overwritten by 3; if it does not, throws a warning
      }

      hold <- hnew
      hnew <- hnew / shrink.factor
      if (hnew < range[1]) hnew <- range[1]
      i <- i + 1
      ulist[[i]] <- getRatio(FUN, x, hnew, ratio.last = ulist[[i-1]], ratio.beforelast = ulist[[i-2]])
    }

    if (any(!ulist[[i]]$monotone)) break
  }
  hopt <- ulist[[i-1]]$h # No monotonicity = bad
  hprev <- ulist[[i-2]]$h
  # Error codes ordered by severity
  if (abs(hopt / hprev - 1) < seq.tol) exitcode <- 2
  if (abs(hopt / range[1] - 1) < seq.tol) {
    exitcode <- 3
    close.left <- TRUE
  }
  if (abs(range[2] / hopt - 1) < seq.tol) exitcode <- 3
  if (i >= maxit) exitcode <- 4

  # !!! If exitcode e, return the last one
  msg <- switch(exitcode + 1,
                "successfully found a monotonicity violation",
                "", # The code cannot be 1 here
                "step size did not change between iterations",
                paste0("step size too close to the ", if (close.left)
                  "left" else "right", " end of the range; consider starting from a ",
                  if (close.left) "larger" else "smaller", " h0 value."),
                "maximum number of iterations reached")

  if (exitcode == 2) warning(paste0("The step size did not change between iterations. ",
                             "This should not happen. Send a bug report to andrei.kostyrka@gmail.com."))
  if (exitcode == 3 && !close.left)
    warning(paste0("The algorithm terminated at the right range of allowed step sizes. ",
                   "Either h0 is too low and the bisection step overshot the next value, ",
                   "or h0 was too large and the truncation error estimate is invalid. ",
                   "Please try a slightly larger and a slightly smaller h0."))

  if (diagnostics) {
    diag.list <- list(h = do.call(c, lapply(ulist, "[[", "h")),
                      x = do.call(rbind, lapply(ulist, "[[", "x")),
                      f = do.call(rbind, lapply(ulist, "[[", "f")),
                      deriv = do.call(c, lapply(ulist, "[[", "deriv")),
                      est.error = do.call(rbind, lapply(ulist, "[[", "est.error")),
                      monotone = do.call(rbind, lapply(ulist, "[[", "monotone")))
  }

  ret <- list(par = hopt, value = ulist[[i-1]]$deriv,
              counts = c(preliminary = i.prelim, main = i - i.prelim),
              exitcode = exitcode, message = msg,
              abs.error = sum(ulist[[i-1]]$est.error),
              iterations = if (diagnostics) diag.list else NULL)

  if (abs(x) + 7e-6 < ret$par)
    warning(paste0("The found step size, ", ret$par, ", exceeds the absolute value of x, ",
                   abs(x), ". Can FUN handle opposite-sign arguments? This seems unreliable. ",
                   "Try a different starting value h0 to be sure."))

  return(ret)
}

#' Automatic step selection for numerical derivatives
#'
#' @param x Numeric vector or scalar: the point at which the derivative is computed
#'   and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric vector or scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param method Character to choose between \insertCite{curtis1974choice}{pnd},
#' modified Curtis--Reid,
#'   \insertCite{dumontet1977determination}{pnd}, and \insertCite{stepleman1979adaptive}{pnd}.
#' @param method.args A named list of tuning parameters for the method. If \code{NULL},
#'   default values are used. See the documentation for the respective methods. Note that
#'   if \code{method.args$diagnostics} is \code{TRUE}, full iteration history
#'   including all function evaluations is returned; different methods have
#'   slightly different diagnostic outputs.
#' @param cores Integer specifying the number of parallel processes to use. Recommended
#'   value: the number of physical cores on the machine minus one.
#' @param ... Passed to FUN.
#'
#' @return A list similar to the one returned by \code{optim()} and made of
#'   concatenated individual elements coordinate-wise lists: \code{par} -- the optimal
#'   step sizes found, \code{value} -- the estimated numerical gradient,
#'   \code{counts} -- the number of iterations for each coordinate,
#'   \code{abs.error} -- an estimate of the total approximation error
#'   (sum of truncation and rounding errors),
#'   \code{exitcode} -- an integer code indicating the termination status:
#'   \code{0} indicates optimal termination within tolerance,
#'   \code{1} means that the truncation error (CR method) or the third derivative
#'   (DV method) is zero and large step size is preferred,
#'   \code{2} is returned if there is no change in step size within tolerance,
#'   \code{3} indicates a solution at the boundary of the allowed value range,
#'   \code{4} signals that the maximum number of iterations was reached.
#'   \code{message} -- summary messages of the exit status.
#'   If \code{method.ards$diagnostics} is \code{TRUE}, \code{iterations} is a list of lists
#'   including the full step size search path, argument grids, function values on
#'   those grids, estimated error ratios, and estimated derivative values for
#'   each coordinate.
#' @export
#'
#' @seealso [step.CR()] for Curtis--Reid (1974) and its modification,
#'   [step.DV()] for Dumontet--Vignes (1977), and
#'   [step.SW()] for Stepleman--Winarsky (1979).
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' gradstep(x = 1, FUN = sin, method = "CR")
#' gradstep(x = 1, FUN = sin, method = "CRm")
#' gradstep(x = 1, FUN = sin, method = "DV")
#' gradstep(x = 1, FUN = sin, method = "SW")
#' # Works for gradients
#' gradstep(x = 1:4, FUN = function(x) sum(sin(x)))
gradstep <- function(x, FUN, h0 = 1e-5 * (abs(x) + (x == 0)),
                     method = c("SW", "CR", "CRm", "DV"),
                     method.args = NULL, cores = 2, ...) {
  # TODO: implement "all"
  method <- method[1]
  dot.args <- list(...)
  if (is.null(h0)) h0 <- 1e-5 * (abs(x) + (x == 0))
  f0 <- .safeF(FUN, x, ...)
  if (length(f0) > 1) stop("The function FUN must return a scalar.")
  if (is.na(f0)) stop(paste0("Could not compute the function value at ", x, ". FUN(x) must be finite."))
  if (length(x) == 1 && length(h0) > 1) stop("The argument 'h0' must be a scalar for scalar 'x'.")
  if (length(x) > 1 && length(h0) == 1) h0 <- rep(h0, length(x))
  if (length(x) != length(h0)) stop("The argument 'h0' must have length 1 or length(x).")
  # The h0 and range arguments are updated later
  default.args <- list(CR = list(h0 = h0[1], version = "original", aim = 100, acc.order = 2, tol = 10,
                                 range = h0[1] / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4, diagnostics = FALSE),
                       CRm = list(h0 = h0[1], version = "modified", aim = 1, acc.order = 2, tol = 4,
                                  range = h0[1] / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4, diagnostics = FALSE),
                       DV = list(h0 = h0[1], range = h0[1] / c(1e6, 1e-6), alpha = 4/3,
                                 ratio.limits = c(1/15, 1/2, 2, 15), maxit = 40L, diagnostics = FALSE),
                       SW = list(h0 = h0[1], shrink.factor = 2, range = h0[1] / c(1e12, 1e-8),
                                 seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2,
                                 maxit = 40L, diagnostics = FALSE))
  margs <- default.args[[method]]
  if (!is.null(method.args)) {
    bad.args <- setdiff(names(method.args), names(margs))
    if (length(bad.args) > 0) {
      stop(paste0("The following arguments are not supported by the ", method, " method: ",
                  .pasteAnd(bad.args)))
    }
    margs[names(method.args)] <- method.args
  }
  conflicting.args <- intersect(names(margs), names(dot.args))
  if (length(conflicting.args) > 0)
    stop(paste0("The arguments ", .pasteAnd(conflicting.args), " of your function coincide with ",
           " the arguments of the ", method, " method. Please write a wrapper for FUN that would ",
           "incorporate the '...' explicitly."))
  autofun <- switch(method, CR = step.CR, CRm = step.CR, DV = step.DV, SW = step.SW)
  if (length(x) == 1) {
    f.args <- c(margs, x = x, FUN = FUN)
    ret <- do.call(autofun, f.args)
  } else {
    f.arg.list <- lapply(seq_along(x), function(i) {
      FUN1 <- function(z) { # Scalar version of FUN
        xx <- x
        xx[i] <- z
        FUN(xx, ...)
      }
      margs1 <- margs
      margs1$h0 <- h0[i]
      margs1$range <- if (!is.null(method.args$range)) method.args$range else
        h0[i] / switch(method, CR = c(1e5, 1e-5), CRm = c(1e5, 1e-5),
                       DV = c(1e6, 1e-6), SW = c(1e12, 1e-8))
      return(c(margs1, x = unname(x[i]), FUN = FUN1))
    })
    ret.list <- lapply(f.arg.list, function(arg1) do.call(autofun, arg1))
    ret <- list(par = do.call(c, lapply(ret.list, "[[", "par")),
                value = do.call(c, lapply(ret.list, "[[", "value")),
                counts = do.call(rbind, lapply(ret.list, "[[", "counts")),
                exitcode = do.call(c, lapply(ret.list, "[[", "exitcode")),
                message = do.call(c, lapply(ret.list, "[[", "message")),
                abs.err = do.call(c, lapply(ret.list, "[[", "abs.error")),
                iterations = lapply(ret.list, "[[", "iterations"))
    ret[names(ret) != "counts"] <- lapply(ret[names(ret) != "counts"],
                                          function(z) {names(z) <- names(x); z})
    if (method == "SW") rownames(ret$counts) <- names(x) else names(ret$counts) <- names(x)
  }
  return(ret)
}
