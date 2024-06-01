#' Gradient and Jacobian computation with parallel capabilities
#'
#' Computes numerical derivatives, gradients and Jacobians. Supports both two-sided
#' (central) and one-sided (forward or backward) derivatives. Calculations can be
#' executed on multiple cores to cut down execution time for slow functions or
#' to attain higher accuracy faster.
#'
#' @param FUN A function returning a numeric scalar or a vector.
#'   If the function returns a vector, the output will be is a Jacobian.
#'   If instead of \code{FUN}, \code{func} is passed, as in \code{numDeriv::grad},
#'   it will be reassigned to \code{FUN} with a warning.
#' @param x Numeric vector or scalar: point at which the derivative is estimated.
#'   \code{FUN(x)} must return a finite value.
#' @param h Numeric scalar or vector specifying the step size for the numerical
#'   difference. Must be length 1 or match \code{length(x)}.
#'   If too large, the slope of the secant poorly estimates the derivative;
#'   if too small, it leads to numerical instability due to the function value rounding.
#' @param deriv.order Integer indicating the derivative order,
#'   \eqn{\mathrm{d}^m / \mathrm{d}x^m}{d^m/dx^m}.
#' @param acc.order Integer specifying the desired accuracy order.
#'   The error typically scales as \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.
#' @param side Integer scalar or vector indicating difference type:
#'   \code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
#'   Using \code{2} (for 'two-sided') triggers a warning and is treated as \code{0}.
#'   with a warning.
#'   Central differences are recommended unless computational cost is prohibitive.
#' @param f0 Optional numeric scalar or vector: if provided and applicable, used
#'   where the stencil contains zero (i.e. \code{FUN(x)} is part of the sum)
#'   to save time. Currently ignored.
#' @param parallel String indicating the parallel execution mode.
#'   This argument is ignored if a valid \code{cluster} is provided,
#'   as the cluster settings take precedence.
#'   On Mac and Linux, use \code{"fork"} to invoke \code{mclapply} or \code{"PSOCK"}
#'   to create a PSOCK cluster dynamically.
#'   On Windows, use \code{"PSOCK"} if you have not created the cluster yet, although
#'   using a dedicated cluster via \code{cluster} is highly recommended to reduce overhead.
#'   \code{"auto"} switches to \code{"PSOCK"} on Windows (with a warning)
#'   and to \code{"fork"} on everything else.
#' @param cores Integer specifying the number of parallel processes to use.
#' @param cluster An optional cluster object. Ensure that necessary objects and packages
#'   are available in the cluster environment with \code{clusterExport()} and \code{clusterEvalQ()}.
#' @param load.balance Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
#'   or enables load balancing with \code{parLapplyLB()}.
#' @param func Deprecated; for \code{numDeriv::grad()} compatibility only.
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @details
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
#' @examples
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
                 h = .Machine$double.eps^(1 / (deriv.order + acc.order)),
                 f0 = NULL,
                 parallel = c("auto", "fork", "PSOCK", "none"), cores = 2,
                 cluster = NULL, load.balance = TRUE, func = NULL,
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
      stop("Pass the function for differencing as the named argument 'FUN' or the first unnamed argument.")
    }
  }

  # 'side', 'deriv.order', 'acc.order', 'h' must align with the length of x
  if (is.null(side)) side <- numeric(n) # NULL --> default central, 0
  if (length(side) == 1) side <- rep(side, n)
  side[!is.finite(side)] <- 0 # NA --> default 'central
  if (any(is2 <- side == 2)) {
    side[is2] <- 0 # Interpreting '2' as two-sided = central
    warning("Interpreting 'side = 2' as '2-sided central differences'; please use side = 0.")
  }
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

  if (length(deriv.order) == 1) deriv.order <- rep(deriv.order, n)
  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)
  # side is already a vector
  # TODO: the part where step is compared to step.CR, step.DV etc.
  if (length(h) == 1) h <- rep(h, n)
  if (length(deriv.order) != n) stop("The argument 'deriv.order' must have length 1 or length(x).")
  if (length(acc.order) != n) stop("The argument 'acc.order' must have length 1 or length(x).")
  if (length(h) != n) stop("The argument 'h' (step size) must have length 1 or length(x).")

  parallel <- parallel[1]
  parallel <- parallel != "none"
  # TODO: proper clusters

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
  fvals <- if (parallel && cores > 1) {
    parallel::mclapply(X = xvals, FUN = function(x) FUN(x, ...), mc.cores = cores, mc.preschedule = !load.balance)
  } else {
    lapply(xvals, function(x) FUN(x, ...))
  }
  # else if (parallel & !use.mclapply & !is.null(cluster)) {
  #if (load.balance) parallel::parLapplyLB(cl = cluster, X = xx, fun = function(x) func(x, ...)) else
  #    parallel::parLapply(cl = cluster, X = xx, fun = function(x) func(x, ...))

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
  return(jac)
}
