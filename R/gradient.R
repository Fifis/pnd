#' Gradient and Jacobian computation with parallel capabilities
#'
#' Computes numerical derivatives, gradients and Jacobians. Supports both two-sided
#' (central) and one-sided (forward or backward) derivatives. Calculations can be
#' executed on multiple cores to cut down execution time for slow functions or
#' to attain higher accuracy faster. Currently for Mac & Linux only because
#' Windows cannot handle \code{parallel::mclapply()}. A \code{parallel::parLapply}
#' version is in development.
#'
#' @param FUN A function returning a numeric scalar or a vector.
#'   If the function returns a vector, the output will be is a Jacobian.
#'   If instead of \code{FUN}, \code{func} is passed, as in \code{numDeriv::grad},
#'   it will be reassigned to \code{FUN} with a warning.
#' @param x Numeric vector or scalar: point at which the derivative is estimated.
#'   \code{FUN(x)} must return a finite value.
#' @param h Numeric scalar, vector, or character specifying the step size for the numerical
#'   difference. If character (\code{"CR"}, \code{"CRm"}, \code{"DV"}, or \code{"SW"}),
#'   calls \code{gradstep()} with the appropriate step-selection method.
#'   Must be length 1 or match \code{length(x)}.
#' @param deriv.order Integer indicating the derivative order,
#'   \eqn{\mathrm{d}^m / \mathrm{d}x^m}{d^m/dx^m}.
#' @param acc.order Integer specifying the desired accuracy order.
#'   The error typically scales as \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.
#' @param side Integer scalar or vector indicating difference type:
#'   \code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
#'   Central differences are recommended unless computational cost is prohibitive.
#' @param f0 Optional numeric scalar or vector: if provided and applicable, used
#'   where the stencil contains zero (i.e. \code{FUN(x)} is part of the sum)
#'   to save time. Currently ignored.
#' @param h0 Numeric scalar of vector: initial step size for automatic search with
#'   \code{gradstep()}.
#' @param method.args A named list of tuning parameters passed to \code{gradstep()}.
#' @param cores Integer specifying the number of parallel processes to use. Recommended
#'   value: the number of physical cores on the machine minus one.
#' @param load.balance Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
#'   or enables load balancing with \code{parLapplyLB()}.
#' @param func Deprecated; for \code{numDeriv::grad()} compatibility only.
#' @param report Integer: if \code{0}, returns a gradient without any attributes; if \code{1},
#'   attaches the step size and its selection method: \code{2} or higher, attaches the full
#'   diagnostic output (overrides \code{diagnostics = FALSE} in \code{method.args}).
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @details
#'
#' If the step size is too large, the slope of the secant poorly estimates the derivative;
#' if it is too small, it leads to numerical instability due to the function value rounding.
#'
#' The optimal step size for one-sided differences typically approaches Mach.eps^(1/2)
#' to balance the Taylor series truncation error with the rounding error due to storing
#' function values with limited precision. For two-sided differences, it is proportional
#' to Mach.eps^(1/3). However, selecting the best step size typically requires knowledge
#' of higher-order derivatives, which may not be readily available. Future releases
#' will allow character arguments to invoke automatic data-driven step-size selection.
#'
#' The use of \code{f0} can reduce computation time similar to the use of \code{f.lower}
#' and \code{f.upper} in \code{uniroot()}.
#'
#' @return Depends on the output of \code{FUN}. If \code{FUN} returns a scalar:
#'   returns a gradient vector matching the length of \code{x}. If \code{FUN} returns a vector:
#'   returns a Jacobian matrix with dimensions \code{length(FUN(x)), length(x)}.
#'   Unlike the output of \code{numDeriv::grad} and \code{numDeriv::jacobian},
#'   this output preserves the names of \code{x} and \code{FUN(x)}.
#'
#' @seealso [gradstep()] for automatic step-size selection.
#'
#' @examples
#'
#' f <- function(x) sum(sin(x))
#' g1 <- Grad(x = 1:4, FUN = f)
#' g2 <- Grad(x = 1:4, FUN = f, h = 7e-6)
#' g2 - g1 # Tiny differences due to different step sizes
#' g3 <- Grad(x = 1:4, FUN = f, h = "SW")
#' g3.full <- Grad(x = 1:4, FUN = f, h = "SW", report = 2)
#' attr(g3.full, "step.search")$exitcode  # Success
#'
#' \dontrun{
#' slowFun <- function(x) {Sys.sleep(0.05); print(x, digits = 12); sum(sin(x))}
#' slowFunVec <- function(x) {Sys.sleep(0.05); print(x, digits = 12)
#'                            c(sin = sum(sin(x)), exp = sum(exp(x)))}
#' true.g <- cos(1:4)  # Analytical gradient
#' true.j <- rbind(cos(1:4), exp(1:4)) # Analytical Jacobian
#' x0 <- c(each = 1, par = 2, is = 3, named = 4)
#'
#' # Compare computation times
#' system.time(g.slow <- numDeriv::grad(slowFun, x = x0) - true.g)
#' system.time(j.slow <- numDeriv::jacobian(slowFunVec, x = x0) - true.j)
#' system.time(g.fast <- Grad(slowFun, x = x0, cores = 4) - true.g)
#' system.time(j.fast <- Grad(slowFunVec, x = x0, cores = 4) - true.j)
#' system.time(j.fast4 <- Grad(slowFunVec, x = x0, acc.order = 4, cores = 4) - true.j)
#'
#' # Compare accuracy
#' rownames(j.slow) <- rep("numDeriv.jac", nrow(j.slow))
#' rownames(j.fast) <- paste0("pnd.jac.order2.", rownames(j.fast))
#' rownames(j.fast4) <- paste0("pnd.jac.order.4", rownames(j.fast4))
#' # Discrepancy
#' print(rbind(numDeriv.grad = g.slow, pnd.Grad = g.fast, j.slow, j.fast, j.fast4), 2)
#' # The order-4 derivative is more accurate for functions
#' # with non-zero third and higher derivatives -- look at pnd.jac.order.4
#'}
#'
#' @export
Grad <- function(FUN, x,
                 deriv.order = 1L, side = 0,
                 acc.order = ifelse(abs(side) == 1, 1L, 2L),
                 h = (abs(x) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order)),
                 h0 = NULL, method.args = list(), f0 = NULL,
                 cores = 2, load.balance = TRUE, func = NULL,
                 report = 1L,
                 ...) {
  n <- length(x)
  if (!(length(side) %in% c(1, n))) stop("The 'side' argument must have length 1 or same length as x.")
  #########################################
  # BEGIN compatibility with numDeriv::grad
  if (missing(FUN)) {
    if (is.function(func)) {
      FUN <- func
      warning("Use the argument 'FUN' to pass the function for differencing to Grad instead of 'func'.")
    } else {
      stop("Pass the function for differencing as the named argument 'FUN' or as the first unnamed argument.")
    }
  }

  # 'side', 'deriv.order', 'acc.order', 'h' must align with the length of x
  if (is.null(side)) side <- numeric(n) # NULL --> default central, 0
  if (length(side) == 1) side <- rep(side, n)
  side[!is.finite(side)] <- 0 # NA --> default 'central
  if (!all(side %in% -1:1)) stop("The 'side' argument must contain values 0 for central, 1 for forward and 2 for backward difference.")

  ell <- list(...)
  if (!is.null(m <- ell$method)) {
    if (identical(m, "simple")) {
      acc.order <- ifelse(side == 0, 2, 1)
    }
    if (identical(m, "Richardson")) {
      if (!is.null(ma <- ell$method.arguments)) {
        if (is.numeric(ma$r)) acc.order <- ma$r
        # TODO: finish grabbing all the arguments
      }
    }
  }
  # END compatibility with numDeriv::grad
  #######################################

  if (!is.function(FUN)) stop("'FUN' must be a function.")
  if (length(deriv.order) == 1) deriv.order <- rep(deriv.order, n)
  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)
  # side is already a vector
  # TODO: the part where step is compared to step.CR, step.DV etc.
  autostep <- FALSE
  if (is.character(h)) {
    method <- h
    if (report == 2) method.args$diagnostics <- TRUE
    h.auto <- gradstep(x = x, FUN = FUN, h0 = h0, method = method, method.args = method.args, ...)
    h <- h.auto$par
    autostep <- TRUE
    # TODO: use this gradient already
  } else if (any(h <= 0)) {
    stop("The argument 'h' (step size) must be positive.")
  }
  if (length(h) == 1) h <- rep(h, n)
  if (length(deriv.order) != n) stop("The argument 'deriv.order' must have length 1 or length(x).")
  if (length(acc.order) != n) stop("The argument 'acc.order' must have length 1 or length(x).")
  if (length(h) != n) stop("The argument 'h' (step size) must have length 1 or length(x).")

  xlist <- lapply(1:n, function(i) {
    sw <- fdCoef(deriv.order = deriv.order[i], acc.order = acc.order[i], side = side[i])
    b <- sw$stencil
    w  <- sw$weights
    bh <- sw$stencil * h[i]
    dx <- matrix(0, ncol = n, nrow = length(b))
    dx[, i] <- bh
    xmat <- matrix(rep(x, length(b)), ncol = n, byrow = TRUE) + dx
    xmat <- as.data.frame(xmat)
    xmat$index <- i
    xmat$weights <- w
    xmat
  })
  xdf <- do.call(rbind, xlist)
  xm <- as.matrix(xdf[, 1:n])
  colnames(xm) <- names(x)
  xvals <- lapply(seq_len(nrow(xdf)), function(i) xm[i, ])

  # Parallelising the task in the most efficient way possible, over all values of all grids
  fvals <- if (cores > 1) {
    parallel::mclapply(X = xvals, FUN = function(x) FUN(x, ...), mc.cores = cores, mc.preschedule = !load.balance)
  } else {
    lapply(xvals, function(x) FUN(x, ...))
  }
  if (any(!sapply(fvals, function(x) is.numeric(x) | is.na(x))))
    stop("'FUN' must output numeric values only, but non-numeric values were returned.")

  # The output can be vector-valued, which is why we ensure that dimensions are not dropped
  fvals <- do.call(rbind, fvals)
  wf <- fvals * xdf$weights
  wf <- split(as.data.frame(wf), f = xdf$index) # Matrices lose dimensions when split
  wf <- lapply(wf, function(x) colSums(as.matrix(x)))
  jac <- unname(do.call(cbind, wf))
  jac <- sweep(jac, 2, h^deriv.order, "/")
  if (nrow(jac) == 1) {
    jac <- drop(jac)
    if (!is.null(names(x))) names(jac) <- names(x)
  } else {
    if (!is.null(names(x))) colnames(jac) <- names(x)
    if (!is.null(colnames(fvals))) rownames(jac) <- colnames(fvals)
  }

  if (report > 0) {
    attr(jac, "step.size") <- h
    if (autostep) {
      attr(jac, "step.size.method") <- method
    } else if (all(h == (abs(x) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order)))) {
      attr(jac, "step.size.method") <- "default"
    } else {
      attr(jac, "step.size.method") <- "user-supplied"
    }
    if (autostep & report > 1) attr(jac, "step.search") <- h.auto
  }

  return(jac)
}
