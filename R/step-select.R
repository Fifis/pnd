#' Default step size at given points
#'
#' Compute an appropriate finite-difference step size based on the magnitude of x,
#' derivation order, and accuracy order. If the function and its higher derivatives
#' belong to the same order of magnitude, this step is near-optimal. For small x,
#' returns a hard bound to prevent large machine-rounding errors.
#'
#' @inheritParams GenD
#'
#' @return A numeric vector of the same length as `x` with positive step sizes.
#' @export
#'
#' @examples
#' stepx(10^(-10:2))
#' stepx(10^(-10:2), deriv.order = 2, acc.order = 4)
stepx <- function(x, deriv.order = 1, acc.order = 2, zero.tol = sqrt(.Machine$double.eps)) {
  x <- abs(x)
  n <- length(x)
  i1 <- x < zero.tol
  i3 <- x > 1
  i2 <- (!i1) & (!i3)
  ret <- rep(zero.tol, length(x))

  if (length(deriv.order) == 1) deriv.order <- rep(deriv.order, n)
  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)
  if (length(deriv.order) != n) stop("The argument 'deriv.order' must have length 1 or length(x).")
  if (length(acc.order) != n) stop("The argument 'acc.order' must have length 1 or length(x).")

  ad <- deriv.order + acc.order
  if (any(i3)) ret[i3] <- x[i3] * .Machine$double.eps^(1/ad[i3])
  if (any(i2)) {  # Exponential interpolation via log(y) ~ a + b*log(x)
    # f(sqrt(macheps)) = macheps^0.5, f(1) = macheps^(1/(a+d))
    leps <- log(.Machine$double.eps)
    a <- leps / ad[i2]
    b <- 1 - a / leps * 2
    ret[i2] <- exp(a + b * log(x[i2]))
  }
  return(ret)
}

# Internal functions to get the quantities used in computing derivatives
# The idea is based on evaluating expressions in parallel with properly exported variables,
# e.g. from the dots (...) argument, as in the optimParallel package.
# Returns a generator object: a list with cluster and a function that can run locally.
getValsCR <- function(FUN, x, h, max.rel.error, vanilla, cores, cl, preschedule, ...) {
  xgrid <- if (vanilla) x + c(-h, 0, h) else x + c(-h, -h/2, h/2, h)
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)

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
  eround <- max(abs(fgrid)) * max.rel.error / h
  ratio <- etrunc / eround
  deriv <- if (vanilla) c(cd = cd, bd = bd, fd = fd) else c(cd = cd, cd.half = cd.half, cd4 = cd4)
  ret <- list(h = h, x = xgrid, f = fgrid, ratio = ratio, deriv = deriv,
              est.error = c(trunc = etrunc, round = eround))
  return(ret)
}

# Generator for Dumontet--Vignes
getValsDV <- function(FUN, x, k, max.rel.error, cores, cl, preschedule, ...) {
  xgrid <- x + c(-2*k, -k, k, 2*k)
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)

  tgrid <- fgrid * c(-0.5, 1, -1, 0.5) # T1, ..., T4 from the paper
  A <- sum(tgrid[tgrid > 0]) # Only positive terms
  B <- sum(tgrid[tgrid < 0]) # Only negative terms
  # TODO: error if fgrid has different sign
  f3sup <- (A / (1-max.rel.error) + B / (1+max.rel.error)) / k^3
  f3inf <- (A / (1+max.rel.error) + B / (1-max.rel.error)) / k^3
  f3 <- sum(tgrid) / k^3 # Estimate of third derivative
  L <- f3sup / f3inf
  L <- max(abs(L), abs(1/L))
  ret <- list(k = k, x = xgrid, f = fgrid, ratio = L,
              deriv = c(f3inf = f3inf, f3 = f3, f3sup = f3sup))
  return(ret)
}

# Generator for plugin
getValsPlugin <- function(FUN, x, h, stage, cores, cl, preschedule, ...) {
  pow <- if (stage == 1) 3 else 1
  s <- fdCoef(deriv.order = pow)
  xgrid <- x + s$stencil * h * if (stage == 1) .Machine$double.eps^(-2/15) else 1

  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)

  f0 <- if (stage == 1) mean(fgrid[2:3]) else mean(fgrid)  # An approximation to f(x)
  cd <- sum(fgrid * s$weights) / h^pow

  list(x = xgrid, f = fgrid, cd = cd, f0 = f0)
}


# Generator for Stepleman--Winarsky
getValsSW <- function(FUN, x, h, do.f0 = FALSE, ratio.last = NULL,
                      ratio.beforelast = NULL, max.rel.error,
                      cores, cl, preschedule, ...) {
  xgrid <- if (do.f0) x + c(-h, 0, h) else x + c(-h, h)

  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)

  cd <- (fgrid[length(fgrid)] - fgrid[1]) / (2*h)

  # Richardson extrapolation to estimate the truncation error
  etrunc <- if (is.null(ratio.last)) NA_real_ else abs((cd - ratio.last$deriv) / (1 - (ratio.last$h / h)^2))
  eround <- max.rel.error * max(abs(fgrid)) / h

  # Check monotonicity of the last three derivatives
  dlast <- if (!is.null(ratio.last) && !is.null(ratio.beforelast)) {
    c(ratio.beforelast$deriv, ratio.last$deriv, cd)
  } else {
    NULL
  }
  # Are derivative changes all positive or all negative?
  monotone.fp <- if (!is.null(dlast)) prod(diff(dlast)) > 0 else NA
  monotone.dfp <- if (!is.null(dlast)) abs(dlast[2] - dlast[3]) <= abs(dlast[2] - dlast[1]) else NA

  list(h = h, x = if (do.f0) xgrid[c(1, 3)] else xgrid, f = if (do.f0) fgrid[c(1, 3)] else fgrid,
       deriv = cd, f0 = if (do.f0) fgrid[2] else NULL, est.error = c(trunc = etrunc, round = eround),
       monotone = c(monotone.fp, monotone.dfp))
}

# Generator for Mathur
getValsM <- function(FUN, x, cores, cl, preschedule, ...) {
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = x, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)
  return(matrix(fgrid, ncol = 2))
}

#' Curtis--Reid automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param max.rel.error Error bound for the relative function-evaluation error
#'   (\eqn{\frac{\hat f(\hat x) - f(x)}{f(x)}}{(^f(^x) - f(x))/f(x)}). Measures how noisy a function is.
#'   If the function is relying on numerical optimisation routines, consider setting to
#'   \code{sqrt(.Machine$double.eps)}.
#'   If the function has full precision to the last bit, set to \code{.Machine$double.eps/2}.
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
#' @inheritParams runParallel
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
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences;
#'       especially useful for computationally expensive functions).
#'     \item \code{counts} – the number of iterations (each iteration includes three function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination within tolerance.
#'         \item \code{1} – Third derivative is zero; large step size preferred.
#'         \item \code{2} – No change in step size within tolerance.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Maximum number of iterations reached.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the full step size search path, argument grids,
#'       function values on those grids, estimated error ratios, and estimated derivative values.
#'   }
#' @export
#'
#' @details
#' The arguments passed to \code{...} must not partially match those of [step.CR()]. For example, if
#' \code{cl} exists, then, attempting to avoid cluster export by using
#' \code{step.CR(f, x, h = 1e-4, cl = cl, a = a)} will result in an error: \code{a} matches \code{aim}
#' and \code{acc.order}. Redefine the function for this argument to have a name that is not equal
#' to the beginning of one of the arguments of [step.CR()].
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.CR(x = 2, f)
#' step.CR(x = 2, f, h0 = 1e-3)
#' step.CR(x = 2, f, version = "modified")
#' step.CR(x = 2, f, version = "modified", acc.order = 4)
#'
#' # A bad start: too far away
#' step.CR(x = 2, f, h0 = 1000)  # Bad exit code + a suggestion to extend the range
#' step.CR(x = 2, f, h0 = 1000, range = c(1e-10, 1e5))  # Problem solved
#'
#' library(parallel)
#' cl <- makePSOCKcluster(names = 2, outfile = "")
#' abc <- 2
#' f <- function(x, abc) {Sys.sleep(0.02); abc*sin(x)}
#' x <- pi/4
#' system.time(step.CR(f, x, h = 1e-4, cores = 1, abc = abc))  # To remove speed-ups
#' system.time(step.CR(f, x, h = 1e-4, cores = 2, abc = abc))  # Faster
#' f2 <- function(x) f(x, abc)
#' clusterExport(cl, c("f2", "f", "abc"))
#' system.time(step.CR(f2, x, h = 1e-4, cl = cl))  # Also fast
#' stopCluster(cl)
step.CR <- function(FUN, x, h0 = 1e-5*max(abs(x), sqrt(.Machine$double.eps)),
                    max.rel.error = .Machine$double.eps^(7/8), version = c("original", "modified"),
                    aim = if (version[1] == "original") 100 else 1,
                    acc.order = c(2L, 4L),
                    tol = if (version[1] == "original") 10 else 4,
                    range = h0 / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                    cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                    cl = NULL, ...) {
  if (length(x) != 1) stop("The step-size selection can handle only univariate inputs. ",
                           "For 'x' longer than 1, use 'gradstep'.")
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  version <- match.arg(version)
  vanilla <- (version == "original")
  acc.order <- acc.order[1]
  cores <- min(cores, if (vanilla) 3 else 4)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)
  if (range[1] < 2*.Machine$double.eps) range[1] <- 2*.Machine$double.eps  # Avoiding zeros
  if (tol <= 1) stop("The tolerance 'tol' must be >=1 (e.g. 4).")
  if (acc.order == 4 && vanilla) {
    warning("The 'original' Curtis--Reid 1974 algorithm does not support 4th-order accuracy. Using acc.order = 2.")
    acc.order <- 2
  }
  if (acc.order == 4) aim <- aim * .Machine$double.eps^(-2/15)
  target <- sort(c(aim / tol, aim * tol))

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  iters <- list()
  i <- 1
  exitcode <- 0

  while (i <= maxit) {
    hold <- if (i > 1) hnew else NA
    hnew <- if (i > 1) hold * (aim / max(iters[[i-1]]$ratio, 1))^(if (vanilla) 1/2 else 1/3) else h0

    # Check the relative change of the step size, which is possible at
    # the 2nd iteration even before the 2nd error calculation
    if (i > 1) {
      if (abs(hnew/hold - 1) < seq.tol) {
        exitcode <- 2
        break  # Step 4: if the step size does not change, stop
      } else { # 4a, b: outside the range, replace with the border
        if (hnew < range[1]) hnew <- range[1]
        if (hnew > range[2]) hnew <- range[2]
        if (max(hnew/hold, hold/hnew) - 1 < seq.tol) {
          exitcode <- 3
          break
        }
      }
    }

    res.i <- getValsCR(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error,
                       vanilla = vanilla, cores = cores, cl = cl, preschedule = preschedule, ...)
    iters[[i]] <- res.i
    if (any(bad <- !is.finite(res.i$f))) {
      bad.iters <- 0
      while (TRUE) {
        bad.iters <- bad.iters + 1
        hnew <- hnew / 2
        if (hnew < max(range[1], .Machine$double.eps))
          stop("step.CR: Could not compute the function value at ", pasteAnd(res.i$x[bad]),
               " after ", bad.iters, " attempts of step shrinkage",
               ".\nChange the range, which is currently [", pasteAnd(range),
               "], and/or\ntry a different starting h0, which is currently ", h0, ".")
        res.i <- getValsCR(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error,
                           vanilla = vanilla, cores = cores, cl = cl, preschedule = preschedule, ...)
        if (!any(bad <- !is.finite(res.i$f))) break
      }
      iters[[i]] <- res.i
    }

    if (res.i$ratio >= target[1] && res.i$ratio <= target[2]) {
      break # Successful termination by matching the range
    }
    if (res.i$ratio < .Machine$double.eps^(4/5)) {
      exitcode <- 1  # Zero truncation error -> only rounding error
      hnew <- hnew * 16 # For linear or quadratic functions, a large h is preferred
      iters[[i+1]] <- getValsCR(FUN = FUN, x = x, h = hnew,  max.rel.error = max.rel.error,
                                vanilla = vanilla, cores = cores, cl = cl, preschedule = preschedule, ...)
      break
    }

    i <- i + 1
  }

  i <- length(iters)
  if (i >= maxit) exitcode <- 4

  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation error is close to zero, large step is favoured",
                "step size did not change between iterations",
                paste0("step size landed on the range ", if (res.i$h == range[1]) "left" else
                         "right", " end; consider extending the range"),
                "maximum number of iterations reached")

  diag.list <- list(h = do.call(c, lapply(iters, "[[", "h")),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    deriv = do.call(rbind, lapply(iters, "[[", "deriv")),
                    est.error = do.call(rbind, lapply(iters, "[[", "est.error")),
                    ratio = do.call(c, lapply(iters, "[[", "ratio")))

  ret <- list(par = res.i$h,
              value = if (acc.order == 4) unname(res.i$deriv["cd4"]) else unname(res.i$deriv["cd"]),
              counts = i, exitcode = exitcode, message = msg,
              abs.error = res.i$est.error,
              method = if (vanilla) "Curtis--Reid" else "Modified Curtis--Reid",
              iterations = diag.list)
  class(ret) <- "stepsize"
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
#' @param max.rel.error Positive numeric scalar > 0 indicating the maximum relative
#'   error of function evaluation. For highly accurate functions with all accurate bits
#'   is equal to \code{.Machine$double.eps/2}. For noisy functions (derivatives, integrals,
#'   output of optimisation routines etc.), it is higher, typically
#'   \code{sqrt(.Machine$double.eps)}. Dumontet and Vignes recommend
#'   \code{.Machine$double.eps^(3/4) = 2e-12} for common functions.
#' @param ratio.limits Numeric vector of length 2 defining the acceptable ranges
#'   for step size: the algorithm stops if the relative perturbation of the third derivative by
#'   amplified rounding errors falls within this range.
#' @param maxit Maximum number of algorithm iterations to avoid infinite loops in cases
#'   the desired relative perturbation factor cannot be achieved within the given \code{range}.
#'   Consider extending the range if this limit is reached.
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{dumontet1977determination}{pnd} algorithm.
#' If the estimated third derivative is exactly zero, the function assumes a third
#' derivative of 1 to prevent division-by-zero errors.
#'
#' Note: the iteration history tracks the third derivative, not the first.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (each iteration includes four function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination within tolerance.
#'         \item \code{1} – Third derivative is zero; large step size preferred.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Maximum number of iterations reached; optimal step size is within the allowed range.
#'         \item \code{5} – Maximum number of iterations reached; optimal step size
#'           was outside allowed range and had to be snapped to a boundary.
#'         \item \code{6} – No search was performed (used when \code{maxit = 1}).
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the full step size search path (note: for the third derivative),
#'       argument grids, function values on those grids, and estimated third derivative values.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.DV(x = 2, f)
#' step.DV(x = 2, f, h0 = 1e-3)
#'
#' # Plug-in estimator with only one evaluation of f'''
#' step.DV(x = 2, f, maxit = 1)
#' step.plugin(x = 2, f)
step.DV <- function(FUN, x, h0 = 1e-5*max(abs(x), sqrt(.Machine$double.eps)),
                    range = h0 / c(1e6, 1e-6), max.rel.error = .Machine$double.eps^(7/8),
                    ratio.limits = c(2, 15), maxit = 40L,
                    cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                    cl = NULL, ...) {
  if (length(x) != 1) stop("The step-size selection can handle only univariate inputs. ",
                           "For 'x' longer than 1, use 'gradstep'.")
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 4)
  k0 <- h0 * .Machine$double.eps^(-2/15)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)
  range3 <- range * .Machine$double.eps^(-2/15)

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  iters <- list()
  i <- 1
  exitcode <- 0
  ndownwards <- 0 # For counting the number of downwards shrinkages

  while (i <= maxit) {
    if (i == 1) k <- k0
    res.i <- getValsDV(FUN = FUN, x = x, k = k, max.rel.error = max.rel.error,
                       cores = cores, cl = cl, preschedule = preschedule, ...)
    iters[[i]] <- res.i
    if (any(bad <- !is.finite(res.i$f))) {
      bad.iters <- 0
      while (TRUE) {
        bad.iters <- bad.iters + 1
        k <- k / 2
        if (k < max(range[1], .Machine$double.eps))
          stop("step.DV: Could not compute the function value at ", pasteAnd(res.i$x[bad]),
               " after ", bad.iters, " attempts of step shrinkage",
               ".\nChange the range, which is currently [", pasteAnd(range),
               "], and/or\ntry a different starting h0, which is currently ", h0, ".")
        res.i <- suppressWarnings(getValsDV(FUN = FUN, x = x, k = k, max.rel.error = max.rel.error,
                                  cores = cores, cl = cl, preschedule = preschedule, ...))
        if (!any(bad <- !is.finite(res.i$f))) break
      }
      iters[[i]] <- res.i
    }

    # Quick rule of thumb: stop after the first iteration
    if (maxit == 1) {
      exitcode <- 6
      break
    }

    # If the estimate of f''' is near-zero, then, the algorithm must stop to avoid division by near-0
    if (abs(res.i$deriv["f3"]) < 8 * .Machine$double.eps) {
      exitcode <- 1
      break
    }

    # If the function sign fluctuates, the step size is too wide
    f.opposite.sign    <- !all(diff(sign(res.i$f)) == c(0, 0, 0))
    # If the numerically contaminated derivative estimates blow up, the step size is too small
    flim.opposite.sign <- !all(diff(sign(res.i$deriv)) == c(0, 0))

    if (res.i$ratio > ratio.limits[2] || flim.opposite.sign) {
      # The rounding error is too high; the step size must be increased
      range3[1] <- k
    } else {
      if (res.i$ratio < ratio.limits[1] || f.opposite.sign) {
        # The rounding error is too small; the step size must be decreased
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

  f3 <- if (exitcode != 1) sum(res.i$f * c(-0.5, 1, -1, 0.5)) / k^3 else 1
  f0 <- mean(res.i$f[2:3]) # Approximately f(x); the error is small for small h
  h <- (1.68 * max.rel.error * abs(f0/f3))^(1/3) # Formula 36 from Dumontet & Vignes (1977)

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
  etrunc <- unname(abs(res.i$deriv["f3"])) * h^2 / 6     # Formula for 'em' from Dumontet (1973), 1.4.2 (p. 37)
  eround <- max.rel.error * max(abs(fgrid)) / h  # Formula for 'ed' ibid.

  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation error is zero, large step is favoured",
                "",
                paste0("step size too close to the ", side,
                       " end of the range [", printE(range[1]), ", ",
                       printE(range[2]), "]; consider extending it"),
                "maximum number of iterations reached",
                paste0("maximum number of iterations reached and step size occured on the ",
                       side, " end of the range [", printE(range[1]), ", ",
                       printE(range[2]), "]; consider expanding it"),
                "only one iteration requested; rough values returned")

  diag.list <- list(k = do.call(c, lapply(iters, "[[", "k")),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    ratio = do.call(c, lapply(iters, "[[", "ratio")),
                    deriv3 = do.call(rbind, lapply(iters, "[[", "deriv")))

  ret <- list(par = h, value = cd, counts = i, exitcode = exitcode,
              message = msg, abs.error = c(trunc = etrunc, round = eround),
              method = "Dumontet--Vignes",
              iterations = diag.list)
  class(ret) <- "stepsize"
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
#' \deqn{\sqrt[3]{1.5 \frac{f'(x)}{f'''(x)} \epsilon_{\mathrm{mach}}}.}{(1.5 mach.eps * f' / f''')^(1/3).}
#' If the estimated third derivative is too small, the function assumes a third
#' derivative of 1 to prevent division-by-zero errors.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (in this case, it is 2).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Termination with checks passed tolerance.
#'         \item \code{1} – Third derivative is exactly zero; large step size preferred.
#'         \item \code{2} – Third derivative is too close to zero; large step size preferred.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the two-step size search path, argument grids,
#'       function values on those grids, and estimated third derivative values.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.plugin(x = 2, f)
#' step.plugin(x = 0, f)  # f''' = 0, setting a large one
step.plugin <- function(FUN, x, h0 = 1e-5*max(abs(x), sqrt(.Machine$double.eps)),
                        max.rel.error = .Machine$double.eps^(7/8), range = h0 / c(1e4, 1e-4),
                        cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                        cl = NULL, ...) {
  # TODO: add zero.tol everywhere
  if (length(x) != 1) stop("The step-size selection can handle only univariate inputs. ",
                           "For 'x' longer than 1, use 'gradstep'.")
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 4)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  iters <- vector("list", 2)
  iters[[1]] <- getValsPlugin(FUN = FUN, x = x, h = h0, stage = 1,
                              cores = cores, cl = cl, preschedule = preschedule, ...)
  if (any(bad <- !is.finite(iters[[1]]$f))) {
    bad.iters <- 0
    while (TRUE) {
      bad.iters <- bad.iters + 1
      h0 <- h0 / 2
      if (h0 < max(range[1], .Machine$double.eps))
        stop("step.plugin: Could not compute the function value at ", pasteAnd(iters[[1]]$x[bad]),
             " after ", bad.iters, " attempts of step shrinkage",
             ".\nChange the range, which is currently [", pasteAnd(range),
             "], and/or\ntry a different starting h0, which is currently ", h0, ".")
      iters[[1]] <- getValsPlugin(FUN = FUN, x = x, h = h0, stage = 1,
                                  cores = cores, cl = cl, preschedule = preschedule, ...)
      if (!any(bad <- !is.finite(iters[[1]]$f))) break
    }
  }


  cd3 <- iters[[1]]$cd
  f0 <- iters[[1]]$f0

  exitcode <- 0
  # If the estimate of f''' is near-zero, the step-size estimate may be too large --
  # only the modified one needs not be saved
  me13 <- max.rel.error^(1/3)
  if (abs(cd3) == 0) {#' Plug-in step selection
    exitcode <- 1
    h <- pmax(me13, abs(x) / 128)
    cd3 <- me13^2 * abs(x)
  } else if (max(abs(f0), me13^2) / abs(cd3) > sqrt(1/max.rel.error)) {
    # The ratio of f' to f''' is too large -- safeguard against large steps
    # small values of f0 are truncated to macheps^(2/3) ~ 4e-11
    cd3 <- sqrt(max.rel.error) * max(abs(f0), me13^2)
    h <- pmax(me13, abs(x) / 256)
    exitcode <- 2
  } else {  # Everything is OK
    h <- abs(1.5 * f0/cd3 * max.rel.error)^(1/3)
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

  iters[[2]] <- getValsPlugin(FUN = FUN, x = x, h = h, stage = 2,
                              cores = cores, cl = cl, preschedule = preschedule, ...)

  msg <- switch(exitcode + 1,
                "successfully computed non-zero f''' and f'",
                "truncation error is zero, large step is favoured",
                "truncation error is near-zero, large step is favoured",
                paste0("step size too close to the ", side,
                       " end of the reasonable range [", printE(range[1]), ", ",
                       printE(range[2]), "]"))

  diag.list <- list(h = c(f3 = h0 * .Machine$double.eps^(-2/15), f1 = h),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    deriv = c(cd3 = cd3, cd = iters[[2]]$cd))

  etrunc <- abs(cd3) / 6 * h^2
  eround <- 0.5 * .Machine$double.eps * max(abs(iters[[2]]$f)) / h
  ret <- list(par = h, value = iters[[2]]$cd, counts = 2, exitcode = exitcode,
              message = msg, abs.error = c(trunc = etrunc, round = eround),
              method = "plug-in",
              iterations = diag.list)
  class(ret) <- "stepsize"
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
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{stepleman1979adaptive}{pnd} algorithm.
#'
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (each iteration includes four function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination within tolerance.
#'         \item \code{2} – No change in step size within tolerance.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Maximum number of iterations reached.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the full step size search path, argument grids,
#'       function values on those grids, estimated derivative values, estimated error values,
#'       and monotonicity check results.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4  # The derivative at 1 is 4
#' step.SW(x = 1, f)
#' step.SW(x = 1, f, h0 = 1e-9) # Starting too low
#' # Starting somewhat high leads to too many preliminary iterations
#' step.SW(x = 1, f, h0 = 10)
#' step.SW(x = 1, f, h0 = 1000) # Starting absurdly high
#'
#' f <- sin  # The derivative at pi/4 is sqrt(2)/2
#' step.SW(x = pi/4, f)
#' step.SW(x = pi/4, f, h0 = 1e-9) # Starting too low
#' step.SW(x = pi/4, f, h0 = 0.1) # Starting slightly high
#' # The following two example fail because the truncation error estimate is invalid
#' step.SW(x = pi/4, f, h0 = 10)   # Warning
#' step.SW(x = pi/4, f, h0 = 1000) # Warning
step.SW <- function(FUN, x, h0 = 1e-5 * (abs(x) + (x == 0)),
                    shrink.factor = 0.5, range = h0 / c(1e12, 1e-8),
                    seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2, maxit = 40L,
                    cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                    cl = NULL, ...) {
  if (length(x) != 1) stop("The step-size selection can handle only univariate inputs. ",
                           "For 'x' longer than 1, use 'gradstep'.")
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 3)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  i <- 1
  exitcode <- 0
  main.loop <- close.left <- FALSE
  first.main <- TRUE
  iters <- list()

  while (i <= maxit) {
    if (!main.loop) {
      if (i == 1) {
        res.i <- getValsSW(FUN = FUN, x = x, h = h0, max.rel.error = max.rel.error, do.f0 = TRUE,
                           ratio.last = NULL, ratio.beforelast = NULL,
                           cores = cores, cl = cl, preschedule = preschedule, ...)
        iters[[i]] <- res.i
        f0 <- res.i$f0  # f(x)
        f <- res.i$f    # f(x +- h)
        hnew <- res.i$h
      }
      if (!is.finite(f0)) stop("Could not compute the function value at ", x, ". FUN(x) must be finite.")
      if (any(bad <- !is.finite(res.i$f))) {
        bad.iters <- 0
        while (TRUE) {
          bad.iters <- bad.iters + 1
          hnew <- hnew / 2
          if (hnew < max(range[1], .Machine$double.eps))
            stop("step.SW: Could not compute the function value at ", pasteAnd(res.i$x[bad]),
                 " after ", bad.iters, " attempts of step shrinkage",
                 ".\nChange the range, which is currently [", pasteAnd(range),
                 "], and/or\ntry a different starting h0, which is currently ", h0, ".")
          res.i <- getValsSW(FUN = FUN, x = x, h = h0, max.rel.error = max.rel.error, do.f0 = TRUE,
                             ratio.last = NULL, ratio.beforelast = NULL,
                             cores = cores, cl = cl, preschedule = preschedule, ...)
          if (!any(bad <- !is.finite(res.i$f))) break
        }
        iters[[i]] <- res.i
      }

      # First check: are the function values of different signs?
      # If yes, h is too large -- jump to the step-shrinking main loop:
      if (prod(sign(res.i$f)) < 0) { # f(x+h) should be close to f(x-h)
        main.loop <- TRUE
        i.prelim <- i
        next
      } else { # The bisection part
        while (!main.loop && i <= maxit) {
          hold <- hnew
          # N: approx. number of inaccurate digits due to rounding at h0 and hopt
          Nh0 <- log10(abs(f0)) - log10(2) - log10(hnew) - log10(abs(res.i$deriv))
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
            res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                               ratio.last = iters[[i-1]], ratio.beforelast = NULL,
                               cores = cores, cl = cl, preschedule = preschedule, ...)
            iters[[i]] <- res.i
            if (abs(hnew/hold - 1) < seq.tol) { # The algorithm is stuck at one range end
              main.loop <- TRUE
              exitcode <- 2
              break
            }

            bad <- !is.finite(res.i$f)
            if (any(bad) && !rounding.small) {  # Not in the original paper, but a necessary fail-safe
              warning("Could not compute the function value at [", pasteAnd(printE(res.i$x[bad])),
                       "]. FUN(", pasteAnd(printE(x)), ") is finite -- try the initial step h0 larger than ",
                       printE(h0), " but smaller than ", printE(hold), ". Halving from ",
                       printE(hnew), " to ", printE(hnew/2), ").")
              for (i in 1:maxit) {
                hnew <- hnew/2
                res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                                   ratio.last = iters[[i-1]], ratio.beforelast = NULL,
                                   cores = cores, cl = cl, preschedule = preschedule, ...)
                iters[[i]] <- res.i
                if (is.finite(res.i$f)) break
              }
              if (!is.finite(res.i$f))
                stop("Could not compute the function value at [", pasteAnd(printE(res.i$x[bad])),
                     "]. FUN(", pasteAnd(printE(x)), ") is finite -- halving did not help.",
                     " Try 'gradstep(..., method = \"K\")' or 'step.K(...)' for a more reliable algorithm.")
            }

            if (any(bad) && !rounding.nottoosmall) {
              warning("Could not compute the function value at ", pasteAnd(res.i$x[bad, , drop = FALSE]),
                      ". FUN(", x, ") is finite -- try a step h0 smaller than ", hnew, ". ",
                      "Halving from ", printE(hnew), " to ", printE(hnew/2), ").")
              for (i in 1:maxit) {
                hnew <- hnew/2
                res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                                   ratio.last = iters[[i-1]], ratio.beforelast = NULL,
                                   cores = cores, cl = cl, preschedule = preschedule, ...)
                iters[[i]] <- res.i
                if (is.finite(res.i$f)) break
              }
              if (!is.finite(res.i$f))
                stop("Could not compute the function value at [",
                     pasteAnd(printE(res.i$x[bad, , drop = FALSE])),
                     "]. FUN(", pasteAnd(printE(x)), ") is finite -- halving did not help.",
                     " Try 'gradstep(..., method = \"M\")' or 'step.M(...)' for a more reliable algorithm.")
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
        res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                           ratio.last = iters[[i-1]], ratio.beforelast = if (j == 1) NULL else iters[[i-2]],
                           cores = cores, cl = cl, preschedule = preschedule, ...)

        iters[[i]] <- res.i
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
      res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                         ratio.last = iters[[i-1]], ratio.beforelast = iters[[i-2]],
                         cores = cores, cl = cl, preschedule = preschedule, ...)
      iters[[i]] <- res.i
    }

    if (!all(res.i$monotone)) break
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

  if (exitcode == 2) warning("The step size did not change between iterations. ",
                             "This should not happen. Send a bug report to https://github.com/Fifis/pnd/issues")
  if (exitcode == 3 && !close.left)
    warning("The algorithm terminated at the right range of allowed step sizes. ",
            "Possible reasons: (1) h0 is too low and the bisection step overshot ",
            "the next value; (2) h0 was too large and the truncation error estimate ",
            "is invalid; (3) the range is too narrow. Please try a slightly larger ",
            "and a slightly smaller h0, or expand the range.")

  if (hopt > 0.01*abs(x) && abs(x) > sqrt(.Machine$double.eps)) {
    exitcode <- 3
    warning("The found step size, ", hopt, ", exceeds 1% of |x|, ",
            abs(x), ", where x is not too small. FUN might poorly behave at x+h ",
            "and x-h due to large steps. Try a different starting value h0 to be sure. ",
            "Returning 0.01|x| .")
    i <- i + 1
    hopt <- 0.01*abs(x)
    res.i <- getValsSW(FUN = FUN, x = x, h = hopt, max.rel.error = max.rel.error,
                       do.f0 = FALSE, ratio.last = if (i > 1) iters[[i-1]] else NULL,
                       ratio.beforelast = if (i > 2) iters[[i-2]] else NULL,
                       cores = cores, cl = cl, preschedule = preschedule, ...)

    iters[[i]] <- res.i
  }

  diag.list <- list(h = do.call(c, lapply(iters, "[[", "h")),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    deriv = do.call(c, lapply(iters, "[[", "deriv")),
                    est.error = do.call(rbind, lapply(iters, "[[", "est.error")),
                    monotone = do.call(rbind, lapply(iters, "[[", "monotone")))


  best.i <- if (exitcode == 3 && !close.left) i-2 else i-1
  ret <- list(par = hopt,
              value = iters[[best.i]]$deriv,
              counts = c(preliminary = i.prelim, main = i - i.prelim),
              exitcode = exitcode, message = msg,
              abs.error = iters[[best.i]]$est.error,
              method = "Stepleman--Winarsky",
              iterations = diag.list)
  class(ret) <- "stepsize"
  return(ret)
}

# The median of the 3 points with the lowest error -- fail-safe return
med3lowest <- function (hgrid, log2etrunc, tstar, hnaive, h0) {
  i.min3 <- rank(log2etrunc, ties.method = "first", na.last = "keep") %in% 1:3
  if (sum(i.min3) >= 3) {
    hopt0 <- sort(hgrid[i.min3])[2]
    i.hopt <- which(hopt0 == hgrid)
    hopt <- hopt0 * (1 / tstar)^(1/3)  # TODO: any power
    exitcode <- 2
  } else {
    hopt0 <- hopt <- hnaive  # At least something should be returned
    i.hopt <- if (sum(i.min3) > 0) min(hgrid[is.finite(log2etrunc)]) else which.min(abs(hgrid-h0))
    exitcode <- 3
  }
  return(list(hopt0 = hopt0, hopt = hopt, i.hopt = i.hopt, exitcode = exitcode))
}

#' Mathur's AutoDX-like automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param max.rel.error Error bound for the relative function-evaluation error
#'   (\eqn{\frac{\hat f(\hat x) - f(x)}{f(x)}}{(^f(^x) - f(x))/f(x)}). Measures how noisy a function is.
#'   If the function is relying on numerical optimisation routines, consider setting to
#'   \code{sqrt(.Machine$double.eps)}.
#'   If the function has full precision to the last bit, set to \code{.Machine$double.eps/2}.
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
#' @inheritParams runParallel
#' @param plot Logical: if \code{TRUE}, plots the estimated truncation and round-off
#'   errors.
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{mathur2012analytical}{pnd} algorithm.
#'
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (each iteration includes two function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination due to a sequence of correct reductions.
#'         \item \code{1} – Reductions are slightly outside the tolerance.
#'         \item \code{2} – Tolerances are significantly violated; an approximate minimum is returned.
#'         \item \code{3} – Not enough finite function values; a rule-of-thumb value is returned.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the step and argument grids,
#'       function values on those grids, estimated derivative values, and estimated error values.
#'   }
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
#' step.M(x = pi/2, f, plot = TRUE)  # Bad case -- TODO a fix
#' step.M(x = pi/4, f, plot = TRUE)
#' step.M(x = pi/4, f, h0 = 1e-9) # Starting low
#' step.M(x = pi/4, f, h0 = 1000) # Starting high
#' # where the truncation error estimate is invalid
step.M <- function(FUN, x, h0 = NULL, max.rel.error = .Machine$double.eps^(7/8), range = NULL,
                   shrink.factor = 0.5, min.valid.slopes = 5L, seq.tol = 0.01,
                   correction = TRUE, plot = FALSE,
                   cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                   cl = NULL, ...) {
  if (length(x) != 1) stop("The step-size selection can handle only univariate inputs. ",
                           "For 'x' longer than 1, use 'gradstep'.")
  cores <- checkCores(cores)
  if (is.null(h0)) { # Setting the initial step to a large enough power of 2
    h0 <- 0.01 * max(abs(x), 1)
    h0 <- 2^round(log2(h0))
  }
  h0 <- unname(h0)
  inv.sf <- 1 / shrink.factor

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  if (is.null(range)) {
    ends <- log(c(2^36, 2^(-24)), base = inv.sf)
    ends <- c(floor(ends[1]), ceiling(ends[2]))
    ends <- inv.sf^ends
    range <- h0 / ends
  }
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)
  # Safety checks for the range size
  hnaive <- (abs(x) * (x!=0) + (x==0)) * .Machine$double.eps^(1/3)
  spans <- c(hnaive/range[1], range[2]/hnaive)
  if (min(spans) < 2^16) {
    range.old <- range
    if (spans[1] < 2^16) range[1] <- hnaive / 2^16
    if (spans[2] < 2^16) range[2] <- hnaive * 2^16
    warning("The initial range [", pasteAnd(printE(range.old)), "] was extended to [",
            pasteAnd(printE(range)), "] to ensure a large-enough search space.")
  }

  exitcode <- 0

  # Creating a sequence of step sizes for evaluation
  sf.sugg <- max(0.5, round(sqrt(shrink.factor), 2))
  range.sugg <- range / c(1024, 1/1024)
  err1 <- paste0("Either increase 'shrink.factor' (e.g. from ", shrink.factor,
                 " to ", sf.sugg, ") to have a finer grid, or increase 'range' (",
                 "from [", pasteAnd(printE(range)), "] to ",
                 "[", pasteAnd(printE(range.sugg)), "]).")
  hgrid <- inv.sf^(floor(log(range[1], base = inv.sf)):ceiling(log(range[2], base = inv.sf)))
  tstar <- (1 + inv.sf) / (1 - shrink.factor^2)  # TODO: n = 2...
  n <- length(hgrid)
  xgrid <- x + c(-hgrid, hgrid)

  fgrid <- getValsM(FUN = FUN, x = xgrid, cores = cores, cl = cl, preschedule = preschedule, ...)

  # TODO: instead of subtracting one, add one
  fplus <- fgrid[, 2]
  fminus <- fgrid[, 1]
  xgrid <- matrix(xgrid, ncol = 2)
  cd <- (fplus - fminus) / hgrid * 0.5
  fd1 <- cd[-1] # Larger step
  fd2 <- cd[-n] # Smaller step
  # Formula 3.7 from Mathur (2012) or (25) from Mathur (2013)
  # Cn*h1^2 = abs((fd2 - fd1) / (1 - (h2/h1)^2), but h2/h1 = const = 1 / shrink.factor
  etrunc <- c(NA, abs((fd2 - fd1) / (1 - shrink.factor^2)))
  log2etrunc <- suppressWarnings(log2(etrunc))
  log2etrunc[!is.finite(log2etrunc)] <- NA
  ldetrunc <- c(NA, diff(log2etrunc))
  signs <- sign(ldetrunc)
  signs.rle <- do.call(data.frame, rle(signs))

  # Rounding error for later
  delta <- .Machine$double.eps / 2 # Error of h|f'(x)true - f'(x)| / |f(x)true|
  fmax <- pmax(abs(fplus), abs(fminus))
  fw <- fdCoef(deriv.order = 1, side = 0, acc.order = 2)  # TODO: generalise later
  f.eps <- max.rel.error * (abs(fw$weights[1]*fminus) + abs(fw$weights[2]*fplus))
  f.delta <- delta*fmax
  eround <- (f.eps + f.delta) / hgrid

  # Find the first sequence of necessary length where the error is increasing
  first.good <- which(signs.rle$values == 1 & signs.rle$lengths > min.valid.slopes)[1]
  if (length(first.good) > 0 && is.finite(first.good)) {
    n.good  <- signs.rle$lengths[first.good]
    i.end <- sum(signs.rle$lengths[1:first.good])
    i.start   <- i.end - n.good + 1
    # Checking the closeness to the slope to the approximation order, e.g. n=2
    i.increasing <- i.start:i.end
    slopes <- ldetrunc[i.increasing]
    valid.h <- hgrid[i.increasing]
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
      i.good <- i.increasing[good.slopes]
    }

    med.slope <- round(stats::median(slopes), 2)
    if (!any(good.slopes)) {
      okay.slopes <- abs(slopes - 2) / 2 <= max(0.1, min(seq.tol * 3, 0.9))
      if (sum(rle(okay.slopes)$values) > 1) okay.slopes <- removeFirstBad(okay.slopes)
      if (any(okay.slopes)) {
        good.slopes <- okay.slopes
        exitcode <- 1
        warning("The estimated truncation error has a slightly wrong reduction rate (~",
                med.slope, ", but should be ~2). ", err1)
      } else {
        warning("The estimated truncation error has a wrong reduction rate (~", med.slope,
                ", but should be ~2). ", err1)
      }
      i.okay <- i.increasing[okay.slopes]
    }
    # Finding the smallest step size with a valid slope and correcting it
    hopt0 <- valid.h[which(good.slopes)[1]]
    if (is.na(hopt0)) {
      m3 <- med3lowest(hgrid, log2etrunc, tstar, hnaive, h0)
      i.hopt <- m3$i.hopt
      hopt0 <- min(h0, m3$hopt0)  # Fail-safe: the bandwidth cannot be too large
      hopt <- min(h0, m3$hopt)
      exitcode <- m3$exitcode
    } else {
      i.hopt <- which(hgrid == hopt0)
      hopt <- hopt0 * (1 / tstar)^(1/3) # TODO: any power
    }
  } else {  # Fail-safe return for cases like x^2 at 0
    m3 <- med3lowest(hgrid, log2etrunc, tstar, hnaive, h0)
    i.hopt <- m3$i.hopt
    hopt0 <- min(h0, m3$hopt0)  # Because hopt0 can be gigantic
    hopt <- min(h0, m3$hopt)
    exitcode <- m3$exitcode

    if (exitcode == 2)
      warning("Could not find a sequence of of ", min.valid.slopes, " reductions ",
              "of the truncation error. Visualise by adding 'plot = TRUE'. ", err1,
              " Finally, try setting 'min.valid.slopes' to 4 or even 3. For now, ",
              "returning the approximate argmin of the total error.")
    if (exitcode == 3)
      warning("There are <3 finite function values on the grid. ",
              "Try setting 'shrink.factor' to 0.9 (close to 1) or checking why the ",
              "function does not return finite values on the range x+[",
              pasteAnd(printE(range)), "] with ", n,
              " exponentially spaced points. Returning a very rough value that may ",
              "not even yield a finite numerical derivative.")
  }

  # !!! If exitcode e, return the last one
  msg <- switch(exitcode + 1,
                "successfully found the optimal step size",
                "successfully found the optimal step size but allowed inaccurate slopes",
                "truncation error reduction rate is too wrong, returning the approximate best step",
                "Fewer than 3 finite function values on the grid, returning the naive step")

  diag.list <- list(h = hgrid, x = xgrid, f = fgrid, deriv = cd,
                    est.error = rbind(NA, cbind(etrunc = etrunc, eround = eround)))

  # TODO: remove the first NA from the output
  ret <- list(par = if (correction) hopt else hopt0, value = cd[i.hopt], counts = n,
              exitcode = exitcode, message = msg,
              abs.error = diag.list$est.error[i.hopt, ],
              method = "Mathur",
              iterations = diag.list)

  if (plot) {
    if (!exists("slopes")) i.increasing <- NULL
    if (!exists("i.good"))
      i.good <- if (!is.null(i.increasing) && exists("good.slopes")) i.increasing[good.slopes] else NULL
    if (!exists("i.okay"))
      i.okay <- if (!is.null(i.increasing) && exists("okay.slopes")) i.increasing[okay.slopes] else NULL
    i.round <- 1:min(i.hopt, i.good, i.okay, i.increasing)
    elabels <- rep("b", n)
    elabels[i.increasing] <- "i"
    elabels[i.round] <- "r"
    elabels[i.good] <- "g"
    elabels[i.okay] <- "o"  # In increasing order of priority

    plotTE(hgrid = hgrid, etotal = etrunc, eround = eround, hopt = hopt,
           elabels = elabels, epsilon = max.rel.error)
  }
  class(ret) <- "stepsize"
  return(ret)
}

#' Kink-based step selection (experimental!)
#'
#' Optimal step-size search using the full range of practical error estimates and
#' numerical optimisation to find the spot where the theoretical total-error shape
#' is best described by the data, and finds the step size where the ratio of
#' rounding-to-truncation error is optimal.
#'
#' @inheritParams step.M
#' @inheritParams GenD
#'
#' @details
#' This function computes the optimal step size for central differences using the statistical
#' kink-search approach.
#' The optimal step size is determined as the minimiser of the total error, which for central
#' differences is V-shaped with the left-branch slope equal to the negative derivation order
#' and the right-branch slope equal to the accuracy order. For standard simple central
#' differences, the slopes are -1 and 2, respectively. The algorithm uses the
#' least-median-of-squares (LMS) penalty and searches for the optimal position of the check
#' that fits the data the best on a bounded 2D rectangle using derivative-free (Nelder--Mead)
#' optimisation.
#' Unlike other algorithms, if the estimated third derivative is too small, the function shape
#' will be different, and two checks are made for the existence of two branches.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (each iteration includes two function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination; a minimum of the V-shaped function was found.
#'         \item \code{1} – Third derivative is too small or noisy; a fail-safe value is returned.
#'         \item \code{2} – Third derivative is nearly zero; a fail-safe value is returned.
#'         \item \code{3} – There is no left branch of the V shape; a fail-safe value is returned.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the step and argument grids,
#'       function values on those grids, estimated derivative values, estimated error values,
#'       and predicted model-based errors.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' step.K(sin, 1, plot = TRUE)
#' step.K(exp, 1, range = c(1e-12, 1e-0), plot = TRUE)
#' step.K(atan, 1, plot = TRUE)
#'
#' # Edge case 1: function symmetric around x0, zero truncation error
#' step.K(sin, pi/2, plot = TRUE)
#' step.K(sin, pi/2, shrink.factor = 0.8, plot = TRUE)
#'
#' # Edge case 1: the truncation error is always zero and f(x0) = 0
#' suppressWarnings(step.K(function(x) x^2, 0, plot = TRUE))
#' # Edge case 2: the truncation error is always zero
#' step.K(function(x) x^2, 1, plot = TRUE)
#' step.K(function(x) x^4, 0, plot = TRUE)
#' step.K(function(x) x^4, 0.1, plot = TRUE)
#' step.K(function(x) x^6 - x^4, 0.1, plot = TRUE)
#' step.K(atan, 3/4, plot = TRUE)
#' step.K(exp, 2, plot = TRUE, range = c(1e-16, 1e+1))
step.K <- function(FUN, x, h0 = NULL, deriv.order = 1, acc.order = 2,
                   range = NULL, shrink.factor = 0.5,
                   max.rel.error = .Machine$double.eps^(7/8), plot = FALSE,
                   cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                   cl = NULL, ...) {
  if (length(x) != 1) stop("The step-size selection can handle only univariate inputs. ",
                           "For 'x' longer than 1, use 'gradstep'.")
  if (!is.numeric(shrink.factor) || shrink.factor <= 0 || shrink.factor >= 1)
    stop("'shrink.factor' must be strictly greater than 0 and less than 1. Recommended: 0.5.")
  if (acc.order %% 2 != 0) stop("'acc.order' must be even for central differences.")

  cores <- checkCores(cores)
  if (is.null(h0)) { # Setting the initial step to a reasonably large power of 2
    h0 <- 0.001 * max(abs(x), 1)
    h0 <- 2^round(log2(h0))
  }
  h0 <- unname(h0)
  inv.sf <- 1 / shrink.factor

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  if (is.null(range)) {
    ends <- log(c(2^36, 2^(-24)), base = inv.sf)
    ends <- c(floor(ends[1]), ceiling(ends[2]))
    ends <- inv.sf^ends
    range <- h0 / ends
  }
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)
  # Safety checks for the range size
  hnaive <- (abs(x) * (x!=0) + (x==0)) * .Machine$double.eps^(1/(deriv.order + acc.order))
  spans <- c(hnaive/range[1], range[2]/hnaive)
  if (min(spans) < 2^16) {
    range.old <- range
    if (spans[1] < 2^16) range[1] <- hnaive / 2^16
    if (spans[2] < 2^16) range[2] <- hnaive * 2^16
    warning("The initial range [", pasteAnd(printE(range.old)), "] was extended to [",
            pasteAnd(printE(range)), "] to ensure a large-enough search space.")
  }

  exitcode <- 0
  hat.check <- NULL

  # Creating a sequence of step sizes for evaluation
  hgrid <- inv.sf^(floor(log(range[1], base = inv.sf)):ceiling(log(range[2], base = inv.sf)))
  n <- length(hgrid)
  xgrid <- x + c(-hgrid, hgrid)

  fgrid <- getValsM(FUN = FUN, x = xgrid, cores = cores, cl = cl, preschedule = preschedule, ...)
  xgrid <- matrix(xgrid, ncol = 2)
  fplus <- fgrid[, 2]
  fminus <- fgrid[, 1]

  # Stencil of powers of 1/shrink.factor because the step is shrunk in this manner
  # Rounding near-integers to integers to handle cases such as sqrt(2)^2 - 2 = 2*MachEps
  s2 <- inv.sf^(0:(acc.order/2 + 1))
  near.int <- (s2/round(s2) - 1) < 16 * .Machine$double.eps
  s2[near.int] <- round(s2[near.int])
  s2 <- c(-rev(s2), s2)
  s1 <- s2[2:(length(s2)-1)]
  s0 <- s1[2:(length(s2)-1)]  # For actual first derivatives
  # Standard + high accuracy for an alternative estimate of the truncation term
  # because it is the next acc.order-th derivative in the Taylor expansion
  stc0  <- fdCoef(deriv.order = deriv.order, acc.order = acc.order,   stencil = s0)
  stc1  <- fdCoef(deriv.order = deriv.order+acc.order, acc.order = acc.order,   stencil = s1)
  stc2  <- fdCoef(deriv.order = deriv.order+acc.order, acc.order = acc.order+2, stencil = s2)

  stc.list <- list(stc0, stc1, stc2)
  cds <- vector("list", 2)
  for (i in 1:3) {
    i.spos <- stc.list[[i]]$stencil > 0  # x+h, x+2h, ...
    i.sneg <- stc.list[[i]]$stencil < 0
    wpos  <-     stc.list[[i]]$weights[i.spos]  # Weights for f(...) starting at h: x+h, x+th, ...
    wneg  <- rev(stc.list[[i]]$weights[i.sneg]) # Weights for f(...) starting at -h: x-h, x-th, ...
    fpos <- sapply(1:sum(i.spos), function(j) c(fplus[j:n], rep(NA, j-1)))
    fneg <- sapply(1:sum(i.sneg), function(j) c(fminus[j:n], rep(NA, j-1)))
    # Estimated required and higher-order derivative times h^acc.order
    cds[[i]] <- (colSums(t(fpos) * wpos) + colSums(t(fneg) * wneg)) / hgrid
  }
  etrunc1 <- cds[[2]] * abs(attr(stc1, "remainder.coef"))  #  f''' * <w, b> / (a+d)!
  etrunc2 <- cds[[3]] * abs(attr(stc2, "remainder.coef"))
  etrunc <- etrunc2
  bad.etrunc <- (!is.finite(etrunc)) | (etrunc == 0)
  etrunc[bad.etrunc] <- etrunc1[bad.etrunc]
  etrunc <- abs(etrunc)
  # plot(abs(cds[[2]]), log = "y"); points(abs(cds[[3]]), col = 2)

  # A value to compare to the threshold for Zero Truncation Error
  good.te <- etrunc != 0 & is.finite(etrunc)
  med.te <- if (any(good.te)) stats::median(etrunc[good.te]) else .Machine$double.eps

  # Approximate rounding error for the plot
  # epsilon is the relative error of evaluation of f(x), e.g. maximum noise
  # delta is the error of |f'(x)true - f'(x)| / |f(x)true|
  # Taking not the exact terms but the first necessary n terms because it does not matter for small h
  # and is too erratic in any case
  stc <- fdCoef(deriv.order = deriv.order, acc.order = acc.order)
  l <- length(stc$stencil) / 2
  fmax <- apply(abs(cbind(fpos[, 1:l], fneg[, 1:l])), 1, max)
  f.eps <- max.rel.error * sum(abs(stc$weights)) * fmax
  f.delta <- 0.5 * .Machine$double.eps * fmax
  eround <- (f.eps + f.delta) / hgrid^deriv.order
  # plot(hgrid, eround, log = "xy")

  # Initial value:
  log2h <- log2(hgrid)
  log2e <- suppressWarnings(log2(etrunc))
  log2e[is.infinite(log2e)] <- NA  # Preserving NaNs just in case
  # plot(log2h, log2e)

  # Check-like function of two parameters: horizontal and vertical location
  checkFit <- function (theta, x, a, m) (theta[2]-m*(x-theta[1]))*(x<theta[1]) + (theta[2]+a*(x-theta[1]))*(x>=theta[1])

  # If the function is symmetric / low-degree polynomial, the truncation error is always zero -- return
  # the safe option; otherwise, carry out a search
  is.finite.etrunc <- is.finite(log2e)
  n.good <- sum(is.finite.etrunc)
  is.good.round <- n.good >= 3
  if (is.good.round) {
    # Identifying regions of rounding, valid truncation, and invalid truncation
    # whilst also safeguarding against noisiness; define the slope at a point as
    # the average right and left slope
    i.etrunc.finite <- which(is.finite.etrunc)
    local.slopes0 <- c(NA, sapply(2:n.good, function(i) {  # Backwards-looking slopes
      ii  <- i.etrunc.finite[i-1:0]
      diff(log2e[ii]) / diff(log2h[ii])
    }))
    local.slopes1 <- (local.slopes0 + c(local.slopes0[-1], NA)) * 0.5
    local.slopes1[c(1, n.good)] <- local.slopes0[c(2, n.good)]
    # Aligning the slopes with the original points
    local.slopes <- rep(NA, n)
    local.slopes[i.etrunc.finite] <- local.slopes1
    # text(log2h[i.etrunc.finite], log2e[i.etrunc.finite], labels = round(local.slopes[i.etrunc.finite], 2), cex = 0.75, col = "#BB0000")

    # Estimated minimum: where the slope changes from negative to positive
    i.trunc.valid <- which(abs((local.slopes - acc.order) / acc.order) < 0.1)
    only.right.branch <- FALSE
    if (length(i.trunc.valid) > 0) {
      zero.trunc <- FALSE
      # Extending it because the local regression used +-1 points at the ends
      i.trunc.valid <- sort(unique(c(i.trunc.valid-1, i.trunc.valid, i.trunc.valid+1)))

      # Edge case: if the slope of etrunc is always acc.order, then, there is no check
      if (0 %in% i.trunc.valid) {
        only.right.branch <- FALSE
        i.trunc.valid <- i.trunc.valid[i.trunc.valid > 0]
        exitcode <- 3
      }

      # If there are more than 1 run, pick the longest one
      # However, if it is shorter than 5 (2 points were added artificially; we want a length-3 run), discard it
      itv.runs <- splitRuns(i.trunc.valid)
      itv.len <- lengths(itv.runs)
      if (length(itv.runs) > 1)
        i.trunc.valid <- itv.runs[[which.max(itv.len)]]
      if (length(i.trunc.valid) >= 5) {
        i.round <- seq(1, (min(i.trunc.valid)-1))
      } else {  # If the truncation error is too erratic, treat everything as a rounding error
        zero.trunc <- med.te < max.rel.error
        i.round <- which(etrunc != 0 & is.finite(etrunc))
      }
    } else {  # Maybe the higher-order derivatives are close to zero compared to the function noise?
      zero.trunc <- med.te < max.rel.error
      # If the truncation error is zero, treat everything as a rounding error
      i.round <- which(etrunc != 0 & is.finite(etrunc))
    }

    if ((!only.right.branch) && (!zero.trunc) && length(i.trunc.valid) >= 3) {
      # It makes sense to equate the two errors when there are both branches with enough points
      # Pseudo-Huber loss with knee radius
      phuber <- function(x, r) r^2 * (sqrt(1 + (x/r)^2) - 1)
      penalty <- function(theta, a, m, subset) {
        ehat <- checkFit(theta = theta, x = log2h, a, m)
        resid <- log2e - ehat
        rmad <- stats::median(abs(resid), na.rm = TRUE)
        if (rmad < max.rel.error) rmad <- mean(abs(resid), na.rm = TRUE)
        mean(phuber(resid, r = rmad)[subset], na.rm = TRUE)
      }
      # Objective function of the two parameters only
      sbst <- sort(unique(c(i.trunc.valid, i.round)))
      fobj <- function (theta) penalty(theta, a = acc.order, m = deriv.order, subset = sbst)

      # Constrained optimisation on a reasonable range
      # Initial value: median of 5 lowest points, but not farther than the the rightmost point of the valid truncation range
      theta0 <- c(stats::median(log2h[rank(log2e, ties.method = "first") <= 5]), min(log2e[sbst], na.rm = TRUE))
      if (theta0[1] > log2h[max(i.trunc.valid)]) theta0 <- c(log2h[min(i.trunc.valid)], log2e[min(i.trunc.valid)])

      # ui <- rbind(diag(2), -diag(2))
      # ci <- c( c(min(log2h), min(log2e, na.rm = TRUE) - 1),
      #          -c(max(log2h), unname(quantile(log2e, 0.75, na.rm = TRUE))))
      # opt <- constrOptim(theta = theta0, f = fobj, grad = NULL, ui = ui, ci = ci, control = list(reltol = 1e-6))
      opt <- stats::optim(par = theta0, fn = fobj, method = "BFGS", control = list(reltol = .Machine$double.eps^(1/3)))
      # Unconstrained optimisation works fine; this is 10x faster than constrOptim

      # Left branch: t2 - m*(x-t1), right branch: t2 + a*(x-t1)
      # Working in logarithms: the difference between the two at hopt must be log(d/a)
      # (a+m) * (x-hopt0) = log(d/a) --> hopt = hopt0 / (d/a)^(1/(d+a))
      # This ensures that etrunc/eround = 1/a
      hopt0   <- if (!zero.trunc) 2^opt$par[1] else 0.001 * max(1, abs(x))
      hopt <- hopt0 * (deriv.order/acc.order)^(1/(deriv.order + acc.order))
      hat.check <- 2^checkFit(opt$par, x = log2h, a = acc.order, m = deriv.order)
    } else {  # Zero truncation --> choosing the smallest error that yields a non-zero total-error estimate
      hopt <- hgrid[min(i.etrunc.finite)]
    }
    log2hopt <- log2(hopt)
    # points(log2h[sbst], log2e[sbst], pch = 16)
    # lines(log2h, log2(hat.check), lty = 2)

    if (length(i.trunc.valid) < 5) {
      warning("Fewer than 3 observations in the truncation branch! Resorting to a safe option.")
      # At the optimum, the rounding error should be approx. |(e^(2/3) f^(2/3) f'''^(1/3))/(2^(2/3) 3^(1/3))|
      # However, if f''' ~ 0, we use only the rounding branch and compare it to the optimal value
      # if f''' ~ 1
      exitcode <- 1
    }
  } else {
    zero.trunc <- TRUE
    exitcode <- 2
  }

  if (exitcode %in% c(1, 2)) {
    complete.rows <- apply(fgrid, 1, function(x) all(is.finite(x)))
    f0 <- abs(max(abs(fgrid[which(complete.rows)[1], ])))
    if (!is.finite(f0)) stop("Could not compute the function value at ", x, ".")
    if (f0 < .Machine$double.eps) f0 <- .Machine$double.eps
    expected.eround <- (.Machine$double.eps^2 * f0^2 / 12)^(1/3)
    # Two cases possible: f(x) = x^2 at x = 0 has eround increasing,
    # which implies that a large fail-safe step can be chosen
    f.er <- is.finite(eround)
    mean.sign <- if (sum(f.er) > 1) mean(sign(diff(eround[f.er]))) else 1
    if (mean.sign > 0.5) {  # Preferring a slightly larger step size
      hopt <- 128*stepx(x = x, deriv.order = deriv.order, acc.order = acc.order,
                        zero.tol = .Machine$double.eps^(1/3))
    } else {
      # Otherwise, for an eround that does not grow convincingly, find the step that approximates
      # the theoretical expected value
      hopt  <- hgrid[which.min(abs(eround - expected.eround))]
    }
    log2hopt <- log2(hopt)
    # TODO: generalise; do the theory!
  } else if (exitcode == 3) {  # If the right branch is growing, choose a smaller step size
    hopt <- hgrid[round(stats::quantile(i.trunc.valid, 0.25))]
    log2hopt <- log2(hopt)
  }

  if (plot) {
    if (any(is.finite.etrunc)) {
      elabels <- rep("b", n)
      elabels[i.round] <- "r"
      if (is.good.round) {  # At least 3 points in the truncation branch
        elabels[i.trunc.valid] <- "g"
        slope.rel.err <- (local.slopes - acc.order) / acc.order
        i.trunc.okay <- which(slope.rel.err > 0 & abs(slope.rel.err) <= 0.5) # Positive slope, +- 50% error
        i.trunc.okay <- setdiff(i.trunc.okay, i.trunc.valid)
        i.trunc.okay <- i.trunc.okay[i.trunc.okay > (max(i.round)-3)]  # To avoid random jumps in the rounding part
        elabels[i.trunc.okay] <- "o"
      } else {
        i.trunc.valid <- i.trunc.okay <- NULL
      }
      i.trunc.increasing <- (c(NA, diff(etrunc)) > 0) | (c(diff(etrunc), NA) > 0)
      i.trunc.increasing <- setdiff(which(i.trunc.increasing), c(i.trunc.okay, i.trunc.valid))
      if (length(i.trunc.increasing) > 0) i.trunc.increasing <- i.trunc.increasing[i.trunc.increasing >= (max(i.round)-2)]
      if (length(i.trunc.increasing) > 0 && length(c(i.trunc.valid, i.trunc.okay)) > 0)
        i.trunc.increasing <- i.trunc.increasing[i.trunc.increasing <= (max(i.trunc.valid, i.trunc.okay)-2)]
      elabels[i.trunc.increasing] <- "i"

      etrunc0 <- etrunc
      etrunc0[etrunc0 == 0] <- NA
      plotTE(hgrid = hgrid, etotal = etrunc0, eround = eround, elabels = elabels,
             hopt = hopt, echeck = hat.check, epsilon = max.rel.error, xlim = range)
    } else {
      warning("Cannot plot: there are no reliable non-zero error estimates.")
    }
  }

  # Approximating f'(hopt) by linear interpolation between 2 closest points
  i.closest <- which(rank(abs(log2(hgrid) - log2hopt), ties.method = "first") <= 2)
  h.closest <- hgrid[i.closest]

  # The grid with the largest stencil remains in memory

  cd.closest <- cds[[1]][i.closest]  # Derivatives at 2 closest points
  if (sum(is.na(cd.closest)) == 1) cd.closest <- rep(cd.closest[is.finite(cd.closest)[1]], 2)
  # If there are no finite values, take the value that is the closest to the centre
  if (sum(is.na(cd.closest)) == 2) {
    i.mid <- which(is.finite(cds[[1]]))
    cd.closest <- rep(hgrid[which.min(abs(i.mid - n/2))], 2)
  }
  cd <- stats::approx(x = h.closest, y = cd.closest, xout = hopt)$y
  if (is.na(cd)) cd <- mean(cd.closest)
  if (!zero.trunc && exists("opt")) {  # If optimisation was carried out
    et <- 2^checkFit(opt$par, x = log2hopt, a = acc.order, m = deriv.order)
  } else {
    et <- med.te
  }
  er <- if (!zero.trunc) et * acc.order else mean(eround[i.closest])
  # TODO: fix this formula, double-check the calculations

  msg <- switch(exitcode + 1,
                "target error ratio reached within tolerance",
                "truncation branch slope is close to zero, relying on the expected rounding error",
                "truncation error is near-zero, relying on the expected rounding error",
                "there is no left branch of the combined error, fitting impossible"

  )

  diag.list <- list(h = hgrid, x = xgrid, f = fgrid, deriv = cbind(f = cds[[1]], fhigher = cds[[2]]),
                    est.error = cbind(etrunc = etrunc, eround = eround))

  ret <- list(par = hopt, value = cd, counts = n, exitcode = exitcode, message = msg,
              abs.error = c(trunc = et, round = er), method = "Kink", iterations = diag.list)
  class(ret) <- "stepsize"
  return(ret)
}

#' Estimated total error plot
#'
#' Visualises the estimated truncation error, rounding error, and total error
#' used in automatic step-size selection for numerical differentiation.
#' The plot follows the approach used in Mathur (2012) and other step-selection methods.
#'
#' @param hgrid Numeric vector: a sequence of step sizes used as the horizontal positions
#'   (usually exponentially spaced).
#' @param etotal Numeric vector: estimated combined error at each step size.
#'   This is typically computed by subtracting a more accurate finite-difference approximation from
#'   a less accurate one.
#' @param eround Numeric vector: estimated rounding error at each step size; usually the best guess
#'   or the upper bound is used.
#' @param echeck Numeric vector: estimated V-shaped check, usually from a fit.
#' @param hopt Numeric scalar (optional): selected optimal step size. If provided,
#'   a vertical line is drawn at this value.
#' @param elabels Character vector of the same length as \code{hgrid} containing the following values:
#'   \code{"r", "g", "o", "i", "b"} for **r**ounding, **g**ood truncation, **o**kay truncation,
#'   **i**ncreasing truncation and **b**ad truncation.
#' @param epsilon Numeric scalar: condition error, i.e. the error bound for the accuracy of the evaluated
#' function; used for labelling rounding error assumptions.
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @returns Nothing (invisible null).
#' @export
#'
#' @examples
#'
#' a <- step.K(sin, 1)
#' hgrid <- a$iterations$h
#' etotal <- a$iterations$est.error[, 1]
#' eround <- a$iterations$est.error[, 2]
#' elabels <- c(rep("r", 32), rep("i", 2), rep("g", 12), rep("b", 15))
#' hopt <- a$par
#' plotTE(hgrid, etotal = 2e-12 * hgrid^2 + 1e-14 / hgrid,
#'        eround = 1e-14 / hgrid, hopt = 0.4, i.increasing = 30:45, i.good = 32:45,
#'        abline.round = c(-46.5, -1))
plotTE <- function(hgrid, etotal, eround, hopt = NULL,
                   elabels = NULL, echeck = NULL,
                   epsilon = .Machine$double.eps^(7/8), ...) {
  cols <- c("#7e1fde", "#328d2d", "#d58726", "#ca203a")
  good.inds <- is.finite(etotal) & (etotal != 0)
  if (!any(good.inds)) {
    warning("The truncation error is either exactly 0 or NA. Nothing to plot.")
    return(invisible(NULL))
  }
  yl <- range(etotal[good.inds])
  plot(hgrid[good.inds], etotal[good.inds],
       log = "xy", bty = "n", ylim = yl,
       ylab = "Estimated abs. error in df/dx", xlab = "Step size",
       main = "Estimated error vs. finite-difference step size", ...)
  graphics::mtext(paste0("dotted line: rounding error assuming max. rel. err. < ", printE(epsilon, 1)), cex = 0.8, line = 0.5)
  i.round      <- which(elabels == "r")
  i.good       <- which(elabels == "g")
  i.ok         <- which(elabels == "o")
  i.increasing <- which(elabels == "i")
  i.bad        <- which(elabels == "b")
  if (length(i.round) > 0)
    graphics::points(hgrid[i.round], etotal[i.round], pch = 16, col = cols[1], cex = 0.9)
  if (length(i.good) > 0)
    graphics::points(hgrid[i.good], etotal[i.good], pch = 16, col = cols[2], cex = 0.9)
  if (length(i.ok) > 0)
    graphics::points(hgrid[i.ok], etotal[i.ok], pch = 16, col = cols[3], cex = 0.9)
  if (length(i.increasing) > 0)
    graphics::points(hgrid[i.increasing], etotal[i.increasing], pch = 16, col = cols[4], cex = 0.9)
  if (length(i.bad) > 0)
    graphics::points(hgrid[i.bad], etotal[i.bad], pch = 4, col = "#00000088", cex = 0.9)

  graphics::lines(hgrid, eround, col = "#000000", lty = 3)

  if (!is.null(echeck)) graphics::lines(hgrid, echeck, col = "#FFFFFFBB", lwd = 3)
  if (!is.null(echeck)) graphics::lines(hgrid, echeck)

  if (!is.null(hopt)) graphics::abline(v = hopt, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Rounding", "Trunc. good", "Trunc. fair", "Trunc. growing", "Invalid"),
                   pch = c(16, 16, 16, 16, 4), col = c(cols, "#000000"),
                   box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(invisible(NULL))
}


#' Automatic step selection for numerical derivatives
#'
#' @param x Numeric vector or scalar: the point at which the derivative is computed
#'   and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric vector or scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param method Character indicating the method: \code{"CR"} for \insertCite{curtis1974choice}{pnd},
#'   \code{"CRm"} for modified Curtis--Reid, "DV" for \insertCite{dumontet1977determination}{pnd},
#'   \code{"SW"} \insertCite{stepleman1979adaptive}{pnd}, \code{"M"} for
#'   \insertCite{mathur2012analytical}{pnd}, \code{"K"} for Kostyrka (2026, exerimental),
#'   and \code{"plugin"} for the single-step plug-in estimator.
#' @param control A named list of tuning parameters for the method. If \code{NULL},
#'   default values are used. See the documentation for the respective methods. Note that
#'   full iteration history including all function evaluations is returned, but
#'   different methods have slightly different diagnostic outputs.
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' We recommend using the Stepleman--Winarsky algorithm because it does not suffer
#' from over-estimation of the truncation error in the Curtis--Reid approach
#' and from sensitivity to near-zero third derivatives in the Dumontet--Vignes
#' approach. It really tries multiple step sizes and handles missing
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
#'   \code{iterations} is a list of lists
#'   including the full step size search path, argument grids, function values on
#'   those grids, estimated error ratios, and estimated derivative values for
#'   each coordinate.
#'
#' @order 1
#' @export
#'
#' @seealso [step.CR()] for Curtis--Reid (1974) and its modification,
#'   [step.plugin()] for the one-step plug-in solution,
#'   [step.DV()] for Dumontet--Vignes (1977),
#'   [step.SW()] for Stepleman--Winarsky (1979),
#'   [step.M()] for Mathur (2012), and
#'   [step.K()] for Kostyrka (2026).
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
gradstep <- function(FUN, x, h0 = NULL,
                     method = c("plugin", "SW", "CR", "CRm", "DV", "M", "K"), control = NULL,
                     cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                     cl = NULL, ...) {
  # TODO: implement "all"
  method <- method[1]
  # TODO: test if Mathur accepts the plot
  if (is.null(h0)) {
    deriv.order <- 1  # TODO: take into account later
    acc.order <- 2
    h0 <- stepx(x, deriv.order = deriv.order, acc.order = acc.order)
    if (method %in% c("M", "K")) h0 <- 2^round(log2(128 * h0))
  }
  h0 <- unname(h0)
  cores <- checkCores(cores)
  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  ell <- list(...)
  if (any(names(ell) == "method.args"))
    stop("'method.args' is an argument to control numDeriv::grad(). ",
         "In pnd::gradstep(), pass the list of step-selection method arguments as 'control'.")
  f0 <- safeF(FUN, x, ...)
  if (length(f0) > 1) stop("Automatic step selection works only when the function FUN returns a scalar.")
  if (is.na(f0)) stop("Could not compute the function value at [", pasteAnd(x), "]. FUN(x) must be finite.")
  if (length(x) == 1 && length(h0) > 1) stop("The argument 'h0' must be a scalar for scalar 'x'.")
  if (length(x) > 1 && length(h0) == 1) h0 <- rep(h0, length(x))
  if (length(x) != length(h0)) stop("The argument 'h0' must have length 1 or length(x).")
  # The h0 and range arguments are updated later
  default.args <- list(plugin = list(h0 = h0[1], max.rel.error = .Machine$double.eps^(7/8),
                                     range = h0[1] / c(1e4, 1e-4),
                                     cores = cores, preschedule = preschedule, cl = cl),
                       CR = list(h0 = h0[1], max.rel.error = .Machine$double.eps^(7/8),
                                 version = "original", aim = 100, acc.order = 2, tol = 10,
                                 range = h0[1] / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                                 cores = cores, preschedule = preschedule, cl = cl),
                       CRm = list(h0 = h0[1], max.rel.error = .Machine$double.eps^(7/8),
                                  version = "modified", aim = 1, acc.order = 2, tol = 4,
                                  range = h0[1] / c(1e5, 1e-5), maxit = 20L, seq.tol = 1e-4,
                                  cores = cores, preschedule = preschedule, cl = cl),
                       DV = list(h0 = h0[1], range = h0[1] / c(1e6, 1e-6), max.rel.error = .Machine$double.eps^(7/8),
                                 ratio.limits = c(2, 15), maxit = 40L,
                                 cores = cores, preschedule = preschedule, cl = cl),
                       SW = list(h0 = h0[1], shrink.factor = 0.5, range = h0[1] / c(1e12, 1e-8),
                                 seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2, maxit = 40L,
                                 cores = cores, preschedule = preschedule, cl = cl),
                       M = list(h0 = h0[1], range = h0[1] / 2^c(36, -24), shrink.factor = 0.5,
                                max.rel.error = .Machine$double.eps^(7/8),
                                min.valid.slopes = 5L, seq.tol = 0.01, correction = TRUE,
                                cores = cores, preschedule = preschedule, cl = cl, plot = FALSE),
                       K = list(h0 = h0[1], range = h0[1] / 2^c(36, -24), shrink.factor = 0.5,
                                max.rel.error = .Machine$double.eps^(7/8),
                                plot = FALSE, cores = cores, preschedule = preschedule, cl = cl))
  margs <- default.args[[method]]
  if (!is.null(control)) {
    bad.args <- setdiff(names(control), names(margs))
    if (length(bad.args) > 0) {
      stop("The following arguments are not supported by the ", method, " method: ",
            pasteAnd(bad.args))
    }
    margs[names(control)] <- control
  }
  conflicting.args <- intersect(names(margs), names(ell))
  if (length(conflicting.args) > 0)
    stop("The arguments ", pasteAnd(conflicting.args), " of your function coincide with ",
         "the arguments of the ", method, " method. Please write a wrapper for FUN that would ",
         "incorporate the '...' explicitly.")
  autofun <- switch(method, plugin = step.plugin, CR = step.CR, CRm = step.CR, DV = step.DV,
                    SW = step.SW, M = step.M, K = step.K)

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
                     DV = c(1e6, 1e-6), SW = c(1e12, 1e-8), M = 2^c(36, -24), K = 2^c(36, -24))
    return(c(margs1, x = unname(x[i]), FUN = FUN1))
  })
  ret.list <- lapply(f.arg.list, function(arg1) do.call(autofun, arg1))
  ret <- list(par = do.call(c, lapply(ret.list, "[[", "par")),
              value = do.call(c, lapply(ret.list, "[[", "value")),
              counts = do.call(rbind, lapply(ret.list, "[[", "counts")),
              exitcode = do.call(c, lapply(ret.list, "[[", "exitcode")),
              message = do.call(c, lapply(ret.list, "[[", "message")),
              abs.error = if (length(x) == 1) unlist(lapply(ret.list, "[[", "abs.error")) else do.call(rbind, lapply(ret.list, "[[", "abs.error")),
              method = method,
              iterations = lapply(ret.list, "[[", "iterations"))
  valid.names <- setdiff(names(ret), c("counts", "abs.error", "method"))
  ret[valid.names] <- lapply(ret[valid.names], function(z) {
    if (is.matrix(z)) rownames(z) <- names(x) else
      if (length(z) == length(x)) names(z) <- names(x) else
        stop("Name error in gradstep (please send a bug report).")
    return(z)
  })
  if (!is.null(dim(ret$counts))) rownames(ret$counts) <- names(x) else names(ret$counts) <- names(x)
  if (!is.null(dim(ret$abs.error))) rownames(ret$abs.error) <- names(x)

  class(ret) <- "gradstep"

  return(ret)
}
