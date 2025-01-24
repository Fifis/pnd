#' Curtis--Reid automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param version Character scalar: \code{"original"} for the original 1974 version by
#'   Curtis and Reid; \code{"modified"} for Kostyrka’s 2025 modification, which adds an
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
#' @inheritParams checkCores
#' @inheritParams runParallel
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
#'
#' # A bad start: too far away
#' step.CR(x = 2, f, h0 = 1000)  # Bad exit code + a suggestion to extend the range
#' step.CR(x = 2, f, h0 = 1000, range = c(1e-10, 1e5))  # Problem solved
step.CR <- function(FUN, x, h0 = 1e-5 * (abs(x) + (x == 0)),
                    version = c("original", "modified"),
                    aim = if (version[1] == "original") 100 else 1,
                    acc.order = c(2L, 4L),
                    tol = if (version[1] == "original") 10 else 4,
                    range = h0 / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                    cores = getOption("pnd.cores"), preschedule = getOption("pnd.preschedule"),
                    cl = NULL, diagnostics = FALSE, ...) {
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  version <- version[1]
  if (!(version %in% c("original", "modified"))) stop("step.CR: 'version' must be either 'original' or 'modified'.")
  vanilla <- version == "original"
  cores <- min(cores, if (vanilla) 3 else 4)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  if (range[2] < range[1]) range <- c(range[2], range[1])
  if (range[1] < 2*.Machine$double.eps) range[1] <- 2*.Machine$double.eps
  if (tol <= 1) stop("The tolerance must be a positive number greater than 1 (e.g. 4).")
  acc.order <- acc.order[1]
  if (acc.order == 4 && vanilla) {
    warning("The original 1974 algorithm does not support 4th-order accuracy. Setting acc.order = 2.")
    acc.order <- 2
  }
  if (acc.order == 4) aim <- aim * .Machine$double.eps^(-2/15)
  target <- sort(c(aim / tol, aim * tol))
  cl <- newCluster(cl = cl, cores = cores)

  getRatio <- function(FUN, x, h, vanilla, ...) {
    xgrid <- x + if (vanilla) c(-h, 0, h) else c(-h, -h/2, h/2, h)
    FUN1 <- function(z) .safeF(FUN, z, ...)
    fgrid <- unlist(runParallel(FUN = FUN1, x = xgrid, cl = cl, preschedule = preschedule))
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

  iters <- list()
  i <- 1
  exitcode <- 0

  while (i <= maxit) {
    hold <- if (i > 1) hnew else NA
    hnew <- if (i > 1) hold * (aim / max(iters[[i-1]]$ratio, 1))^(1/3) else h0

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

    iters[[i]] <- getRatio(FUN = FUN, x = x, h = hnew, vanilla = vanilla, ...)
    if (any(bad <- !is.finite(iters[[i]]$f))) {
      stop(paste0("Could not compute the function value at ", .pasteAnd(iters[[i]]$x[bad]),
                  ". Change the range, which is currently [", .pasteAnd(range),
                  "], and/or try a different starting h0, which is currently ", h0, "."))
    }

    if (iters[[i]]$ratio >= target[1] && iters[[i]]$ratio <= target[2]) {
      break # Successful termination by matching the range
    }
    if (iters[[i]]$ratio == 0) { # Zero truncation error -> only rounding error
      exitcode <- 1
      hnew <- hnew * 4 # For linear or quadratic functions, a large h is preferred
      iters[[i+1]] <- getRatio(FUN = FUN, x = x, h = hnew, vanilla = vanilla, ...)
      break
    }

    i <- i + 1
  }

  i <- length(iters)
  if (i >= maxit) exitcode <- 4

  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation error is exactly zero, large step is favoured",
                "step size did not change between iterations",
                paste0("step size landed on the range ", if (iters[[i]]$h == range[1]) "left" else
                         "right", " end; consider extending the range"),
                "maximum number of iterations reached")
  if (diagnostics) {
    diag.list <- list(h = do.call(c, lapply(iters, "[[", "h")),
                      x = do.call(rbind, lapply(iters, "[[", "x")),
                      f = do.call(rbind, lapply(iters, "[[", "f")),
                      deriv = do.call(rbind, lapply(iters, "[[", "deriv")),
                      est.error = do.call(rbind, lapply(iters, "[[", "est.error")),
                      ratio = do.call(c, lapply(iters, "[[", "ratio")))
  }
  ret <- list(par = iters[[i]]$h,
              value = if (acc.order == 4) unname(iters[[i]]$deriv["cd4"]) else unname(iters[[i]]$deriv["cd"]),
              counts = i, exitcode = exitcode, message = msg,
              abs.error = sum(iters[[i]]$est.error),
              iterations = if (diagnostics) diag.list else NULL)
  return(ret)
}


#' Dumontet--Vignes automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}). This step
#'   size for first derivarives is internallt translated into the initial step size for third
#'   derivatives by multiplying it by the machine epsilon raised to the power -2/15.
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
#' @inheritParams checkCores
#' @inheritParams runParallel
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
#'   \code{4} signals that the maximum number of iterations was reached and the
#'   found optimal step size belongs to the allowed range,
#'   \code{5} occurs when the maximum number of iterations was reached and the
#'   found optimal step size did belong to the allowed range and had to be snapped
#'   to one end.
#'   \code{6} is used when \code{maxit = 1} and no search was performed.
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
#'
#' # Plug-in estimator with only one evaluation of f'''
#' step.DV(x = 2, f, maxit = 1)
step.DV <- function(FUN, x, h0 = 1e-5 * (abs(x) + (x == 0)),
                    range = h0 / c(1e6, 1e-6), alpha = 4/3,
                    ratio.limits = c(1/15, 1/2, 2, 15), maxit = 40L,
                    cores = getOption("pnd.cores"), preschedule = getOption("pnd.preschedule"),
                    cl = NULL,  diagnostics = FALSE, ...) {
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 4)
  k0 <- h0 * .Machine$double.eps^(-2/15)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  if (range[2] < range[1]) range <- c(range[2], range[1])
  range3 <- range * .Machine$double.eps^(-2/15)
  P <- 2^(log2(.Machine$double.eps/2) / alpha)
  cl <- newCluster(cl = cl, cores = cores)

  getRatio <- function(FUN, x, k, P, ...) {
    xgrid <- x + c(-2*k, -k, k, 2*k)
    FUN1 <- function(z) .safeF(FUN, z, ...)
    fgrid <- unlist(runParallel(FUN = FUN1, x = xgrid, cl = cl, preschedule = preschedule))
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

  iters <- list()
  i <- 1
  exitcode <- 0
  ndownwards <- 0 # For counting the number of downwards shrinkages

  while (i <= maxit) {
    if (i == 1) k <- k0
    iters[[i]] <- getRatio(FUN = FUN, x = x, k = k, P = P, ...)
    if (any(bad <- !is.finite(iters[[i]]$f))) {
      stop(paste0("Could not compute the function value at ", .pasteAnd(iters[[i]]$x[bad]),
                  ". Change the range, which is currently [", .pasteAnd(range),
                  "], and/or try a different starting h0, which is currently ", h0, "."))
    }

    # Quick rule of thumb: stop after the first iteration
    if (maxit == 1) {
      exitcode <- 6
      break
    }

    # If the estimate of f''' is near-zero, then, f_inf and f_sup will have opposite signs
    # The algorithm must therefore stop
    if (abs(iters[[i]]$deriv["f3"]) < 8 * .Machine$double.eps) {
      exitcode <- 1
      break
    }

    # TODO: find an improvement for the ratios of opposite signs

    if (iters[[i]]$ratio < ratio.limits[1] || iters[[i]]$ratio > ratio.limits[4]) {
      # The numerical error is too high
      range3[1] <- k
    } else {
      if (iters[[i]]$ratio > ratio.limits[2] && iters[[i]]$ratio < ratio.limits[3]) {
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
  i <- length(iters)

  f3 <- if (exitcode != 1) sum(iters[[i]]$f * c(-0.5, 1, -1, 0.5)) / k^3 else 1
  f0 <- mean(iters[[i]]$f[2:3]) # Approximately f(x); the error is small for small h
  h <- (1.67 * P * abs(f0/f3))^(1/3) # Formula 36 from Dumontet & Vignes (1977)

  if (h < range[1]) {
    h <- range[1]
    exitcode <- 3
    side <- "left"
  }
  if (h > range[2]) {
    h <- range[2]
    exitcode <- 3
    side <- "right"
  }
  # Was the procedure systematically unsuccsessful?
  if (i >= maxit && maxit > 1) {  # Did it waste many iterations in vain?
    exitcode <- if (h == range[1] || h == range[2]) 5 else 4
    side <- if (ndownwards >= maxit/2) "right" else "left"
  }

  h <- h + x # Minimising the representation error
  h <- h - x
  xgrid <- c(x-h, x+h)
  fgrid <- vapply(xgrid, FUN, ..., FUN.VALUE = numeric(1))
  cd <- (fgrid[2] - fgrid[1]) / h / 2
  err <- P*abs(f0)/h/3 + (exitcode != 1) * (h^5*f3^2/36/P/abs(f0) - h^8*abs(f3)^3/648/P^2/f0^2)

  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation error is zero, large step is favoured",
                "",
                paste0("step size too close to the ", side,
                       " end of the range [", .printE(range[1]), ", ",
                       .printE(range[2]), "]; consider extending it"),
                "maximum number of iterations reached",
                paste0("maximum number of iterations reached and step size occured on the ",
                       side, " end of the range [", .printE(range[1]), ", ",
                       .printE(range[2]), "]; consider expanding it"),
                "only one iteration requested; rough values returned")

  if (diagnostics) {
    diag.list <- list(k = do.call(c, lapply(iters, "[[", "k")),
                      x = do.call(rbind, lapply(iters, "[[", "x")),
                      f = do.call(rbind, lapply(iters, "[[", "f")),
                      ratio = do.call(c, lapply(iters, "[[", "ratio")),
                      deriv3 = do.call(rbind, lapply(iters, "[[", "deriv")))
  }

  ret <- list(par = h, value = cd, counts = i, exitcode = exitcode,
              message = msg, abs.error = err,
              iterations = if (diagnostics) diag.list else NULL)
  return(ret)

}

#' Plug-in step selection
#'
#' @inheritParams step.DV
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' plug-in approach.
#' The optimal step size is determined as the minimiser of the total error, which for central
#' finite differences is (assuming minimal bounds for relative rounding errors)
#' \deqn{\sqrt[3]{1.5 \frac{f'(x)}{f'''(x) \epsilon_{\mathrm{mach}}}}}{[(1.5 mach.eps * f' / f''')^(1/3)]}
#' If the estimated third derivative is too small, the function assumes a third
#' derivative of 1 to prevent division-by-zero errors.
#'
#' @return A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
#'   step size found, \code{value} -- the estimated numerical first derivative (central
#'   differences), \code{counts} -- the number of iterations (here, it is 2),
#'   \code{abs.error} -- an estimate of the total approximation error (sum of truncation and
#'   rounding errors),
#'   \code{exitcode} -- an integer code indicating the termination status:
#'   \code{0} indicates termination with checks passed tolerance,
#'   \code{1} means that the third derivative is exactly zero (large step size preferred),
#'   \code{2} signals that the third derivative is too close to zero (large step size preferred),
#'   \code{3} indicates a solution at the boundary of the allowed value range.
#'   \code{message} is a summary message of the exit status.
#'   If \code{diagnostics} is \code{TRUE}, \code{iterations} is a list
#'   including the two-step size search path, argument grids, function values on those grids,
#'   and estimated 3rd derivative values.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.plugin(x = 2, f)
#' step.plugin(x = 0, f, diagnostics = TRUE)  # f''' = 0, setting a large one
step.plugin <- function(FUN, x, h0 = 1e-5 * (abs(x) + (x == 0)), range = h0 / c(1e4, 1e-4),
                        cores = getOption("pnd.cores"), preschedule = getOption("pnd.preschedule"),
                        cl = NULL,  diagnostics = FALSE, ...) {
  # TODO: add zero.tol everywhere
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 4)
  k0 <- h0 * .Machine$double.eps^(-2/15)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  if (range[2] < range[1]) range <- c(range[2], range[1])
  cl <- newCluster(cl = cl, cores = cores)

  s3 <- fdCoef(deriv.order = 3)
  xgrid3 <- x + s3$stencil*k0
  FUN1 <- function(z) .safeF(FUN, z, ...)
  fgrid3 <- unlist(runParallel(FUN = FUN1, x = xgrid3, cl = cl, preschedule = preschedule))
  f0 <- mean(fgrid3[2:3])
  cd3 <- sum(fgrid3 * s3$weights) / k0^3

  iters <- vector("list", 2)
  iters[[1]] <- list(h = c(f = h0, f3 = k0), x = xgrid3, f = fgrid3, deriv = c(f1 = NA, f3 = cd3))

  exitcode <- 0
  # If the estimate of f''' is near-zero, the step-size estimate may be too large --
  # only the modified one needs not be saved
  me13 <- .Machine$double.eps^(1/3)  # 6e-06 in IEEE754
  if (abs(cd3) == 0) {
    exitcode <- 1
    h <- pmax(me13, abs(x) / 128)
    cd3 <- me13^2 * abs(x)
  } else if (max(abs(f0), me13^2) / abs(cd3) > sqrt(1/.Machine$double.eps)) {
    # The ratio of f' to f''' is too large -- safeguard against large steps
    # small values of f0 are truncated to macheps^(2/3) ~ 4e-11
    cd3 <- sqrt(.Machine$double.eps) * max(abs(f0), me13^2)
    h <- pmax(me13, abs(x) / 256)
    exitcode <- 2
  } else {  # Everything is OK
    h <- abs(1.5 * f0/cd3 * .Machine$double.eps)^(1/3)
  }

  if (h < range[1]) {
    h <- range[1]
    exitcode <- 3
    side <- "left"
  } else if (h > range[2]) {
    h <- range[2]
    exitcode <- 3
    side <- "right"
  }

  s1 <- fdCoef(deriv.order = 1)
  xgrid <- x + s1$stencil*h
  fgrid <- unlist(runParallel(FUN = FUN1, x = xgrid, cl = cl, preschedule = preschedule))
  cd <- sum(fgrid * s1$weights) / h

  msg <- switch(exitcode + 1,
                "successfully computed non-zero f''' and f'",
                "truncation error is zero, large step is favoured",
                "truncation error is near-zero, large step is favoured",
                paste0("step size too close to the ", side,
                       " end of the reasonable range [", .printE(range[1]), ", ",
                       .printE(range[2]), "]"))

  iters[[2]] <- list(h = c(f1 = h, f3 = NA), x = xgrid, f = fgrid,
                     deriv = c(f1 = cd, f3 = NA))

  if (diagnostics) {
    diag.list <- list(h = c(f1 = h, f3 = k0),
                      x = do.call(rbind, lapply(iters, "[[", "x")),
                      f = do.call(rbind, lapply(iters, "[[", "f")),
                      deriv = c(cd = cd, cd3 = cd3))
  }

  etrunc <- (mean(fgrid)^2 * abs(cd3) * .Machine$double.eps^2 / 96)^(1/3)
  eround <- max(2*etrunc, (abs(mean(fgrid)) * .Machine$double.eps)^(2/3))
  ret <- list(par = h, value = cd, counts = 2, exitcode = exitcode,
              message = msg, abs.error = c(trunc = etrunc, round = eround),
              iterations = if (diagnostics) diag.list else NULL)
  return(ret)

}


#' Stepleman--Winarsky automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param shrink.factor A scalar less than 1 that is used to multiply the step size
#'   during the search. The authors recommend 0.25, but this may be result in earlier
#'   termination at slightly sub-optimal steps. Change to 0.5 for a more thorough search.
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
#' @inheritParams checkCores
#' @inheritParams runParallel
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
#' # Starting somewhat high leads to too many preliminary iterations
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
step.SW <- function(FUN, x, h0 = 1e-5 * (abs(x) + (x == 0)),
                    shrink.factor = 0.5, range = h0 / c(1e12, 1e-8),
                    seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2, maxit = 40L,
                    cores = getOption("pnd.cores"), preschedule = getOption("pnd.preschedule"),
                    cl = NULL, diagnostics = FALSE, ...) {
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 3)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  if (range[2] < range[1]) range <- c(range[2], range[1])
  cl <- newCluster(cl = cl, cores = cores)

  getRatio <- function(FUN, x, h, do.f0 = FALSE, ratio.last = NULL,
                       ratio.beforelast = NULL, ...) {
    xgrid <- x + if (do.f0) c(-h, 0, h) else c(-h, h)
    FUN1 <- function(z) .safeF(FUN, z, ...)
    fgrid <- unlist(runParallel(FUN = FUN1, x = xgrid, cl = cl, preschedule = preschedule))
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

  i <- 1
  exitcode <- 0
  main.loop <- close.left <- FALSE
  first.main <- TRUE
  iters <- list()

  while (i <= maxit) {
    if (!main.loop) {
      if (i == 1) {
        iters[[i]] <- getRatio(FUN, x, h0, do.f0 = TRUE, ratio.last = NULL, ...)
        f0 <- iters[[i]]$f0
        hnew <- iters[[i]]$h
      }
      if (!is.finite(f0)) {
        stop(paste0("Could not compute the function value at ", x, ". FUN(x) must be finite."))
      }
      if (any(bad <- !is.finite(iters[[i]]$f))) {
        stop(paste0("Could not compute the function value at ", .pasteAnd(iters[[i]]$x[bad]),
                    ". FUN(", x, ") is finite -- reduce the step h0, which is currently ", h0, "."))
      }

      # First check: are the function values of different signs?
      # If yes, h is too large -- jump to the step-shrinking main loop:
      if (prod(sign(iters[[i]]$f)) < 0) { # f(x+h) should be close to f(x-h)
        main.loop <- TRUE
        i.prelim <- i
        next
      } else { # The bisection part
        while (!main.loop && i <= maxit) {
          hold <- hnew
          # N: approx. number of inaccurate digits due to rounding at h0 and hopt
          Nh0 <- log10(abs(f0)) - log10(2) - log10(hnew) - log10(abs(iters[[i]]$deriv))
          Nhopt <- -log10(max.rel.error) / 3
          rounding.nottoosmall <- Nh0 > 0 # Formula 3.16: some digits were lost,
          # the rounding error not zero, the step size not extremely large
          rounding.small <- Nh0 <= Nhopt + log10(shrink.factor) # Formula 3.15:
          # the rounding error is not worse than the truncation error
          if (rounding.small && rounding.nottoosmall) {
            main.loop <- TRUE # We are in the truncation branch -- go to shrinkage directly
            break
          } else { # Either rounding.small (too small) or rounding.nottoosmall (too large)
            i <- i + 1
            # This is what is necessary for 3.15 to hold with a small safety margin
            hnew <- hnew * 10^(Nh0 - (Nhopt + log10(shrink.factor)) + log10(2))
            # If 3.15 does not hold, this increases the step size
            # If 3.16 does not hold, shrink hnew towards the 3.15-optimal value
            if (!rounding.nottoosmall) hnew <- (hold + hnew) / 2 # Bisection

            if (hnew > range[2]) hnew <- range[2] # Safeguarding against huge steps
            iters[[i]] <- getRatio(FUN, x, hnew, do.f0 = FALSE, ratio.last = iters[[i-1]], ...)
            if (abs(hnew/hold - 1) < seq.tol) { # The algorithm is stuck at one range end
              main.loop <- TRUE
              exitcode <- 2
              break
            }

            bad <- !is.finite(iters[[i]]$f)
            if (any(bad) && !rounding.small) {  # Not in the original paper, but a necessary fail-safe
              warning(paste0("Could not compute the function value at [", .pasteAnd(.printE(iters[[i]]$x[bad])),
                             "]. FUN(", .pasteAnd(.printE(x)), ") is finite -- try the initial step h0 larger than ",
                             .printE(h0), " but smaller than ", .printE(hold), ". Halving from ",
                             .printE(hnew), " to ", .printE(hnew/2), ")."))
              for (i in 1:maxit) {
                hnew <- hnew/2
                iters[[i]] <- getRatio(FUN, x, hnew, do.f0 = FALSE, ratio.last = iters[[i-1]], ...)
                if (is.finite(iters[[i]]$f)) break
              }
              if (!is.finite(iters[[i]]$f))
                stop(paste0("Could not compute the function value at [", .pasteAnd(.printE(iters[[i]]$x[bad])),
                             "]. FUN(", .pasteAnd(.printE(x)), ") is finite -- halving did not help.",
                             " Try 'gradstep(..., method = \"M\")' or 'step.M(...)' for a more reliable algorithm."))
            }

            if (any(bad) && !rounding.nottoosmall) {
              warning(paste0("Could not compute the function value at ", .pasteAnd(iters[[i]]$x[bad, , drop = FALSE]),
                             ". FUN(", x, ") is finite -- try a step h0 smaller than ", hnew, ". ",
                             "Halving from ", .printE(hnew), " to ", .printE(hnew/2), ")."))
              for (i in 1:maxit) {
                hnew <- hnew/2
                iters[[i]] <- getRatio(FUN, x, hnew, do.f0 = FALSE, ratio.last = iters[[i-1]], ...)
                if (is.finite(iters[[i]]$f)) break
              }
              if (!is.finite(iters[[i]]$f))
                stop(paste0("Could not compute the function value at [",
                            .pasteAnd(.printE(iters[[i]]$x[bad, , drop = FALSE])),
                             "]. FUN(", .pasteAnd(.printE(x)), ") is finite -- halving did not help.",
                             " Try 'gradstep(..., method = \"M\")' or 'step.M(...)' for a more reliable algorithm."))
            }
          }
        } # End initial step search
      }
      i.prelim <- i # Passing to the first iteration of the main loop with 2 function values saved
    } # End preliminary loop


    if (first.main) { # First main loop: extra h1 and h2 needed
      for (j in 1:2) { # Try a decreasing sequence
        i <- i + 1
        hnew <- hnew * shrink.factor
        if (hnew < range[1]) hnew <- range[1]
        iters[[i]] <- getRatio(FUN, x, hnew, ratio.last = iters[[i-1]],
                               ratio.beforelast = if (j == 1) NULL else iters[[i-2]])
      }
      first.main <- FALSE
    } else { # Monotonicity satisfied, continuing shrinking, only h[i+1] needed
      if (abs(iters[[i]]$h/iters[[i-1]]$h - 1) < seq.tol) {
        exitcode <- 2 # If h did not shrink, it must have hit the lower bound
        break  # or something else went wrong; this code will most likely
        # be overwritten by 3; if it does not, throws a warning
      }

      hold <- hnew
      hnew <- hnew * shrink.factor
      if (hnew < range[1]) hnew <- range[1]
      i <- i + 1
      iters[[i]] <- getRatio(FUN, x, hnew, ratio.last = iters[[i-1]], ratio.beforelast = iters[[i-2]])
    }

    if (any(!iters[[i]]$monotone)) break
  }
  hopt <- iters[[i-1]]$h # No monotonicity = bad
  hprev <- iters[[i-2]]$h
  # Error codes ordered by severity
  if (abs(hopt / hprev - 1) < seq.tol) exitcode <- 2

  # If hprev hits the upper limit, the total error needs to be compared
  if (abs(hprev / range[2] - 1) < seq.tol) {
    abs.error.prev <- unname(iters[[i-1]]$est.error["trunc"] / shrink.factor + iters[[i-2]]$est.error["round"])
    abs.error.opt  <- unname(iters[[i-1]]$est.error["trunc"] + iters[[i-1]]$est.error["round"])
    if (abs.error.prev < abs.error.opt) { # The two-times reduction was unnecessary
      exitcode <- 3
      hopt <- hprev
      iters[[i-2]]$est.error["trunc"] <- unname(iters[[i-1]]$est.error["trunc"] / shrink.factor)
    }
  }
  if (abs(hopt / range[1] - 1) < seq.tol) {
    exitcode <- 3
    close.left <- TRUE
  }

  if (i >= maxit) exitcode <- 4

  # !!! If exitcode e, return the last one
  msg <- switch(exitcode + 1,
                "successfully found a monotonicity violation",
                "", # The code cannot be 1 here
                "step size did not change between iterations",
                paste0("step size too close to the ", if (close.left)
                  "left" else "right", " end of the range; consider extending the range ",
                  "or starting from a ", if (close.left) "larger" else "smaller", " h0 value."),
                "maximum number of iterations reached")

  if (exitcode == 2) warning(paste0("The step size did not change between iterations. ",
                             "This should not happen. Send a bug report to https://github.com/Fifis/pnd/issues"))
  if (exitcode == 3 && !close.left)
    warning(paste0("The algorithm terminated at the right range of allowed step sizes. ",
                   "Possible reasons: (1) h0 is too low and the bisection step overshot ",
                   "the next value; (2) h0 was too large and the truncation error estimate ",
                   "is invalid; (3) the range is too narrow. Please try a slightly larger ",
                   "and a slightly smaller h0, or expand the range."))

  if (diagnostics) {
    diag.list <- list(h = do.call(c, lapply(iters, "[[", "h")),
                      x = do.call(rbind, lapply(iters, "[[", "x")),
                      f = do.call(rbind, lapply(iters, "[[", "f")),
                      deriv = do.call(c, lapply(iters, "[[", "deriv")),
                      est.error = do.call(rbind, lapply(iters, "[[", "est.error")),
                      monotone = do.call(rbind, lapply(iters, "[[", "monotone")))
  }

  best.i <- if (exitcode == 3 && !close.left) i-2 else i-1
  ret <- list(par = hopt,
              value = iters[[best.i]]$deriv,
              counts = c(preliminary = i.prelim, main = i - i.prelim),
              exitcode = exitcode, message = msg,
              abs.error = sum(iters[[best.i]]$est.error),
              iterations = if (diagnostics) diag.list else NULL)

  if (abs(x) + 7e-6 < ret$par)
    warning(paste0("The found step size, ", ret$par, ", exceeds the absolute value of x, ",
                   abs(x), ". It seems unreliable because FUN might poorly behave handle at ",
                   "opposite-sign arguments due to large steps. Try a different starting value h0 to be sure, ",
                   "or ."))

  return(ret)
}

#' Mathur's AutoDX-like automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param shrink.factor A scalar less than 1 that is used to create a sequence of
#'   step sizes. The recommended value is 0.5. Change to 0.25 for a faster search. This
#'   number should be a negative power of 2 for the most accurate representation.
#' @param seq.tol Numeric scalar: maximum relative difference between old and new
#'   step sizes for declaring convergence.
#' @param min.valid.slopes Positive integer: how many points must form a sequence
#'   with the correct slope with relative difference from 2 less than \code{seq.tol}.
#'   If \code{shrink.factor} is small (< 0.33), consider reducing this to 4.
#' @param correction Logical: if \code{TRUE}, returns the corrected step size (last
#'   point in the sequence times a less-than-1 number to account for the possible
#'   continuation of the downwards slope of the total error); otherwise, returns
#'   the grid point that is is lowest in the increasing sequence of valid error
#'   estimates.
#' @inheritParams checkCores
#' @inheritParams runParallel
#' @param diagnostics Logical: if \code{TRUE}, returns the full iteration history
#'   including all function evaluations.
#' @param plot Logical: if \code{TRUE}, plots the estimated truncation and round-off
#'   errors.
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{mathur2012analytical}{pnd} algorithm.
#'
#'
#' @return A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
#'   step size found, \code{value} -- the estimated numerical first derivative (central
#'   differences), \code{counts} -- the number of iterations (each iteration includes
#'   two function evaluations), \code{abs.error} -- an estimate of the total
#'   approximation error (sum of truncation and rounding errors),
#'   \code{exitcode} -- an integer code indicating the termination status:
#'   \code{0} indicates optimal termination due to a sequence of correct reductions,
#'   \code{1} indicates that the reductions are slightly not within tolerance,
#'   \code{2} indicates that the tolerances are so wrong, an approximate minimum is returned,
#'   \code{3} signals that there are not enough finite function values and the rule of thumb is returned.
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
#' step.M(x = 1, f, plot = TRUE)
#' step.M(x = 1, f, h0 = 1e-9) # Starting low
#' step.M(x = 1, f, h0 = 1000) # Starting high
#'
#' f <- sin  # The derivative at pi/4 is sqrt(2)/2
#' step.M(x = pi/2, f)  # Bad case -- TODO a fix
#' step.M(x = pi/4, f, h0 = 1e-9) # Starting low
#' step.M(x = pi/4, f, h0 = 1000) # Starting high
#' # where the truncation error estimate is invalid
step.M <- function(FUN, x, h0 = NULL, range = NULL, shrink.factor = 0.5,
                   min.valid.slopes = 5L, seq.tol = 0.01,
                   correction = TRUE, diagnostics = FALSE, plot = FALSE,
                   cores = getOption("pnd.cores"), preschedule = getOption("pnd.preschedule"),
                   cl = NULL, ...) {
  cores <- checkCores(cores)
  if (is.null(h0)) { # Setting the initial step to a large enough power of 2
    h0 <- 0.01 * (abs(x) + (x == 0))
    h0 <- 2^round(log2(h0))
  }
  h0 <- unname(h0)
  inv.sf <- 1 / shrink.factor
  cl <- newCluster(cl = cl, cores = cores)

  if (is.null(range)) {
    ends <- log(c(2^36, 2^(-24)), base = inv.sf)
    ends <- c(floor(ends[1]), ceiling(ends[2]))
    ends <- inv.sf^ends
    range <- h0 / ends
  }
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  if (range[2] < range[1]) range <- c(range[2], range[1])
  # Safety checks for the range size
  hnaive <- (abs(x) * (x!=0) + (x==0)) * .Machine$double.eps^(1/3)
  spans <- c(hnaive/range[1], range[2]/hnaive)
  if (min(spans) < 2^16) {
    range.old <- range
    if (spans[1] < 2^16) range[1] <- hnaive / 2^16
    if (spans[2] < 2^16) range[2] <- hnaive * 2^16
    warning("The initial range [", .pasteAnd(.printE(range.old)), "] was extended to [",
            .pasteAnd(.printE(range)), "] to ensure a large-enough search space.")
  }

  hgrid <- inv.sf^(floor(log(range[1], base = inv.sf)):ceiling(log(range[2], base = inv.sf)))
  tstar <- (1 + inv.sf) / (1 - shrink.factor^2)  # TODO: n = 2...
  exitcode <- 0
  sf.sugg <- max(0.5, round(sqrt(shrink.factor), 2))
  range.sugg <- range / c(1024, 1/1024)
  err1 <- paste0("Either increase 'shrink.factor' (e.g. from ", shrink.factor,
                 " to ", sf.sugg, ") to have a finer grid, or increase 'range' (",
                 "from [", .pasteAnd(.printE(range)), "] to ",
                 "[", .pasteAnd(.printE(range.sugg)), "]).")

  n <- length(hgrid)
  xgrid <- x + c(-hgrid, hgrid)
  FUN1 <- function(z) .safeF(FUN, z, ...)
  fgrid <- unlist(runParallel(FUN = FUN1, x = xgrid, cl = cl, preschedule = preschedule))
  # TODO: instead of subtracting one, add one
  fplus <- fgrid[(n+1):(2*n)]
  fminus <- fgrid[1:n]
  cd <- (fplus - fminus) / hgrid * 0.5
  fd1 <- cd[-1] # Larger step
  fd2 <- cd[-n] # Small#' @inheritParams runParalleler step
  # Formula 3.7 from Mathur (2012) or (25) from Mathur (2013)
  # Cn*h1^2 = abs((fd2 - fd1) / (1 - (h2/h1)^2), but h2/h1 = const = 1 / shrink.factor
  etrunc <- c(NA, abs((fd2 - fd1) / (1 - shrink.factor^2)))
  log2etrunc <- suppressWarnings(log2(etrunc))
  log2etrunc[!is.finite(log2etrunc)] <- NA
  ldetrunc <- c(NA, diff(log2etrunc))
  signs <- sign(ldetrunc)
  signs.rle <- do.call(data.frame, rle(signs))

  # Rounding error for later
  eps <- .Machine$double.eps / 2
  delta <- .Machine$double.eps / 2 # Error of h|f'(x)true - f'(x)| / |f(x)true|
  fmax <- pmax(abs(fplus), abs(fminus))
  fw <- fdCoef(deriv.order = 1, side = 0, acc.order = 2)  # TODO: generalise later
  f.eps <- eps * (abs(fw$weights[1]*fminus) + abs(fw$weights[2]*fplus))
  f.delta <- delta*fmax
  eround <- (f.eps + f.delta) / hgrid

  # Find the first sequence of necessary length where the error is increasing
  first.good <- which(signs.rle$values == 1 & signs.rle$lengths > min.valid.slopes)[1]
  if (is.finite(first.good) && length(first.good) > 0) {
    n.good  <- signs.rle$lengths[first.good]
    i.end <- sum(signs.rle$lengths[1:first.good])
    i.start   <- i.end - n.good + 1
    # Checking the closeness to the slope to the approximation order, e.g. n=2
    slopes <- ldetrunc[i.start:i.end]
    valid.h <- hgrid[i.start:i.end]
    good.slopes <- abs(slopes - 2) / 2 <= seq.tol  # TODO: generalise with (d)
    # TODO: debug this function, test with shrink.factor = 0.25; the slope seems to be 4

    removeFirstBad <- function(slopes) {
      slopes.rle <- rle(slopes)
      i.last.good.slope <- max(which(slopes))
      i.last.rle.true <- max(which(slopes.rle$values))
      i.first.good.slope <- sum(slopes.rle$lengths[1:(i.last.rle.true-1)]) + 1
      slopes[-(i.first.good.slope:i.last.good.slope)] <- FALSE
      return(slopes)
    }
    if (sum(rle(good.slopes)$values) > 1) {  # More that 1 TRUE sequence, e.g. T F T
      # Mark everything but the last sequence as F
      good.slopes <- removeFirstBad(good.slopes)
    }

    med.slope <- round(stats::median(slopes), 2)
    if (!any(good.slopes)) {
      okay.slopes <- abs(slopes - 2) / 2 <= max(0.1, min(seq.tol * 3, 0.9))
      okay.slopes <- removeFirstBad(okay.slopes)
      if (any(okay.slopes)) {
        good.slopes <- okay.slopes
        exitcode <- 1
        warning(paste0("The estimated truncation error has a slightly wrong reduction rate (~",
                       med.slope, ", but should be ~2). ", err1))
      } else {
        stop(paste0("The estimated truncation error has a wrong reduction rate (~", med.slope,
                    ", but should be ~2). ", err1))
      }
    }
    # Finding the smallest step size with a valid slope and correcting it
    hopt0 <- valid.h[which(good.slopes)[1]]
    i.hopt <- which(hopt0 == hgrid)
    hopt <- hopt0 * (1 / tstar)^(1/3)  # TODO: any power
  } else {
    i.min5 <- rank(log2etrunc, ties.method = "first", na.last = "keep") %in% 1:3
    if (sum(i.min5) >= 3) {
      exitcode <- 2
      warning(paste0("Could not find a sequence of of ", min.valid.slopes, " reductions ",
                     "of the truncation error. Visualise by adding 'plot = TRUE'. ", err1,
                     " Finally, try setting 'min.valid.slopes' to 4 or even 3. For now, ",
                     "returning the approximate argmin of the total error."))
      hopt0 <- sort(hgrid[i.min5])[2]  # The median of the 3 points with the lowest error
      i.hopt <- which(hopt0 == hgrid)
      hopt <- hopt0 * (1 / tstar)^(1/3)  # TODO: any power
    } else {
      warning(paste0("There are ", sum(i.min5), " finite function values on the grid. ",
                     "Try setting 'shrink.factor' to 0.9 (close to 1) or checking why the ",
                     "function does not return finite values on the range x+[",
                     .pasteAnd(.printE(range)), "] with ", n,
                     " exponentially spaced points. Returning a very rough value that may ",
                     "not even yield a numerical derivative."))
      exitcode <- 3
      hopt0 <- hopt <- hnaive  # At least something should be returned
      i.hopt <- if (sum(i.min5) > 0) min(hgrid[is.finite(log2etrunc)]) else which.min(abs(hgrid-h0))
    }
  }

  # !!! If exitcode e, return the last one
  msg <- switch(exitcode + 1,
                "successfully found the optimal step size",
                "successfully found the optimal step size but allowed inaccurate slopes",
                "truncation error reduction rate is too wrong, returning the approximate best step",
                paste0("only ", sum(i.min5), "finite function values on the grid, returning the naive step"))

  diag.list <- list(h = hgrid, x = matrix(xgrid, ncol = 2), f = matrix(fgrid, ncol = 2),
                    deriv = cd,
                    est.error = rbind(NA, cbind(etrunc = etrunc, eround = eround)))

  # TODO: remove the first NA from the output
  ret <- list(par = if (correction) hopt else hopt0, value = cd[i.hopt], counts = n,
              exitcode = exitcode, message = msg,
              abs.error = diag.list$est.error[i.hopt, ],
              iterations = if (diagnostics) diag.list else NULL)

  if (plot) {
    etotal <- etrunc + eround
    evec <- c(etotal, etrunc, eround)
    yl <- range(evec[is.finite(evec) & evec != 0])
    plot(hgrid[etotal != 0], etotal[etotal != 0], log = "xy", bty = "n",
         ylim = yl,
         ylab = "Estimated abs. error in df/dx", xlab = "Step size",
         main = "Estimated error vs. finite-difference step size")
    graphics::points(hgrid[etrunc != 0], etrunc[etrunc != 0], pch = 2, cex = 1.2)
    graphics::points(hgrid[eround != 0], eround[eround != 0], pch = 3)
    graphics::legend("top", c("Truncation", "Rounding", "Total"), pch = c(2, 3, 1), bty = "n", ncol = 3)
    graphics::mtext(paste0("assuming rel. condition err. < ", .printE(eps),
                           ", rel. subtractive err. < ", .printE(delta)), cex = 0.8, line = 0.5)
    if (exists("good.slopes") && any(good.slopes)) {
      good.h <- valid.h[which(good.slopes)]
      i.good.h <- which(hgrid %in% good.h)
      graphics::points(hgrid[i.start:i.end], etrunc[i.start:i.end], lwd = 1.5, col = "#3355FF", cex = 0.80)
      ## TODO: colour okay slopes differently, warn...
    }
    if (exists("i.start")) graphics::points(hgrid[i.start:i.end], etrunc[i.start:i.end],
                                            pch = 16, col = "#CC4422", cex = 0.65)
  }

  return(ret)
}


#' Automatic step selection for numerical derivatives
#'
#' @param x Numeric vector or scalar: the point at which the derivative is computed
#'   and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric vector or scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param zero.tol Small positive integer: if \code{abs(x) >= zero.tol}, then, the automatically
#'   guessed step size is relative (\code{x} multiplied by the step), unless an auto-selection
#'   procedure is requested; otherwise, it is absolute.
#' @param method Character indicating the method: \code{"CR"} for \insertCite{curtis1974choice}{pnd},
#'   \code{"CR"} for modified Curtis--Reid, "DV" for \insertCite{dumontet1977determination}{pnd},
#'   \code{"SW"} \insertCite{stepleman1979adaptive}{pnd}, and "M" for
#'   \insertCite{mathur2012analytical}{pnd}.
#' @param diagnostics Logical: if \code{TRUE}, returns the full iteration history
#'   including all function evaluations. Passed to the appropriate \code{step.XX} function.
#' @param control A named list of tuning parameters for the method. If \code{NULL},
#'   default values are used. See the documentation for the respective methods. Note that
#'   if \code{control$diagnostics} is \code{TRUE}, full iteration history
#'   including all function evaluations is returned; different methods have
#'   slightly different diagnostic outputs.
#' @inheritParams checkCores
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' We recommend using the Mathur algorithm because it does not suffer
#' from over-estimation of the truncation error in the Curtis--Reid approach
#' and from sensitivity to near-zero third derivatives in the Dumontet--Vignes
#' approach. It really tries muliple step sizes simultaneously and handles missing
#' values due to bad evaluations for inadequate step sizes really in a robust manner.
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
#'   [step.DV()] for Dumontet--Vignes (1977),
#'   [step.SW()] for Stepleman--Winarsky (1979), and
#'   [step.M()] for Mathur (2012).
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' gradstep(x = 1, FUN = sin, method = "CR")
#' gradstep(x = 1, FUN = sin, method = "CRm")
#' gradstep(x = 1, FUN = sin, method = "DV")
#' gradstep(x = 1, FUN = sin, method = "SW")
#' gradstep(x = 1, FUN = sin, method = "M")
#' # Works for gradients
#' gradstep(x = 1:4, FUN = function(x) sum(sin(x)))
gradstep <- function(FUN, x, h0 = NULL, zero.tol = sqrt(.Machine$double.eps),
                     method = c("plugin", "SW", "CR", "CRm", "DV", "M"), diagnostics = FALSE, control = NULL,
                     cores = getOption("pnd.cores"), preschedule = getOption("pnd.preschedule"),
                     cl = NULL, ...) {
  # TODO: implement "all"
  method <- method[1]
  # TODO: test if Mathur accepts the plot
  if (is.null(h0)) {
    lttol <- abs(x) < zero.tol
    deriv.order <- 1
    acc.order <- 2
    if (method != "M") {
      h0 <- (abs(x)*(!lttol) + lttol) * 1e-5
    } else {
      h0 <- 2^round(log2(0.01 * (abs(x)*(!lttol) + lttol)))
    }
  }
  cores <- checkCores(cores)
  ell <- list(...)
  if (any(names(ell) == "method.args"))
    stop(paste0("'method.args' is an argument to control numDeriv::grad(). ",
                "In pnd::gradstep(), pass the list of step-selection method arguments as 'control'."))
  f0 <- .safeF(FUN, x, ...)
  if (length(f0) > 1) stop("Automatic step selection works only when the function FUN returns a scalar.")
  if (is.na(f0)) stop(paste0("Could not compute the function value at ", x, ". FUN(x) must be finite."))
  if (length(x) == 1 && length(h0) > 1) stop("The argument 'h0' must be a scalar for scalar 'x'.")
  if (length(x) > 1 && length(h0) == 1) h0 <- rep(h0, length(x))
  if (length(x) != length(h0)) stop("The argument 'h0' must have length 1 or length(x).")
  # The h0 and range arguments are updated later
  default.args <- list(plugin = list(h0 = h0[1], range = h0[1] / c(1e4, 1e-4),
                                 cores = cores, preschedule = preschedule,
                                 cl = cl, diagnostics = diagnostics),
                       CR = list(h0 = h0[1], version = "original", aim = 100, acc.order = 2, tol = 10,
                                 range = h0[1] / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                                 cores = cores, preschedule = preschedule, cl = cl, diagnostics = diagnostics),
                       CRm = list(h0 = h0[1], version = "modified", aim = 1, acc.order = 2, tol = 4,
                                  range = h0[1] / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                                  cores = cores, preschedule = preschedule, cl = cl, diagnostics = diagnostics),
                       DV = list(h0 = h0[1], range = h0[1] / c(1e6, 1e-6), alpha = 4/3,
                                 ratio.limits = c(1/15, 1/2, 2, 15), maxit = 40L,
                                 cores = cores, preschedule = preschedule, cl = cl, diagnostics = diagnostics),
                       SW = list(h0 = h0[1], shrink.factor = 0.5, range = h0[1] / c(1e12, 1e-8),
                                 seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2, maxit = 40L,
                                 cores = cores, preschedule = preschedule, cl = cl, diagnostics = diagnostics),
                       M = list(h0 = h0[1], range = h0[1] / 2^c(36, -24), shrink.factor = 0.5,
                                min.valid.slopes = 5L, seq.tol = 0.01, correction = TRUE,
                                cores = cores, preschedule = preschedule, cl = cl,
                                diagnostics = diagnostics, plot = FALSE))
  margs <- default.args[[method]]
  if (!is.null(control)) {
    bad.args <- setdiff(names(control), names(margs))
    if (length(bad.args) > 0) {
      stop(paste0("The following arguments are not supported by the ", method, " method: ",
                  .pasteAnd(bad.args)))
    }
    margs[names(control)] <- control
  }
  conflicting.args <- intersect(names(margs), names(ell))
  if (length(conflicting.args) > 0)
    stop(paste0("The arguments ", .pasteAnd(conflicting.args), " of your function coincide with ",
           "the arguments of the ", method, " method. Please write a wrapper for FUN that would ",
           "incorporate the '...' explicitly."))
  autofun <- switch(method, plugin = step.plugin, CR = step.CR, CRm = step.CR, DV = step.DV, SW = step.SW, M = step.M)
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
      margs1$range <- if (!is.null(control$range)) control$range else
        h0[i] / switch(method, plugin = c(1e4, 1e-4), CR = c(1e5, 1e-5), CRm = c(1e5, 1e-5),
                       DV = c(1e6, 1e-6), SW = c(1e12, 1e-8), M = 2^c(36, -24))
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
                                          function(z) structure(z, names = names(x)))
    if (method == "SW") rownames(ret$counts) <- names(x) else names(ret$counts) <- names(x)
  }
  return(ret)
}
