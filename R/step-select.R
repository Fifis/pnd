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


# h corresponding to the median of 3 lowest points with a fall-back
med3lowest <- function (y, hgrid, hnaive, tstar) {
  # Hnaive is necessary for index calculation
  i.min3 <- rank(y, ties.method = "first", na.last = "keep") %in% 1:3
  if (sum(i.min3) >= 3) {
    hopt0 <- sort(hgrid[i.min3])[2]
    i.hopt <- which(hgrid == hopt0)
    hopt <- hopt0 * (1 / tstar)^(1/3)  # TODO: any power
    exitcode <- 2L
  } else {
    hopt0 <- hopt <- hnaive  # At least something should be returned
    i.hopt <- if (sum(i.min3) > 0) min(hgrid[is.finite(y)]) else which.min(abs(hgrid - hnaive))
    exitcode <- 3L
  }
  return(list(hopt0 = hopt0, hopt = hopt, i.hopt = i.hopt, exitcode = exitcode))
}


#' Automatic step selection for gradients
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
  if (is.na(f0)) stop("Could not compute the function value at [", toString(x), "]. FUN(x) must be finite.")
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
                                cores = cores, preschedule = preschedule, cl = cl),
                       K = list(h0 = h0[1], range = h0[1] / 2^c(36, -24), shrink.factor = 0.5,
                                max.rel.error = .Machine$double.eps^(7/8),
                                cores = cores, preschedule = preschedule, cl = cl))
  margs <- default.args[[method]]
  if (!is.null(control)) {
    bad.args <- setdiff(names(control), names(margs))
    if (length(bad.args) > 0) {
      stop("The following arguments are not supported by the ", method, " method: ",
            toString(bad.args))
    }
    margs[names(control)] <- control
  }
  conflicting.args <- intersect(names(margs), names(ell))
  if (length(conflicting.args) > 0)
    stop("The arguments ", toString(conflicting.args), " of your function coincide with ",
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
