#' Numerical cross-derivatives with parallel capabilities
#'
#' Computes the derivative of a function with respect to two different arguments.
#' Arbitrary accuracies and sides for different coordinates
#' of the argument vector are supported.
#'
#' @param FUN A function returning a numeric scalar.
#'   If the function returns a vector, the output will be is a Jacobian.
#'   If instead of \code{FUN}, \code{func} is passed, as in \code{numDeriv::grad},
#'   it will be reassigned to \code{FUN} with a warning.
#' @param x Numeric vector or scalar: point at which the derivative is estimated.
#'   \code{FUN(x)} must return a finite value.
#' @param h Numeric scalar, vector, or character specifying the step size for the numerical
#'   difference. If character (\code{"CR"}, \code{"CRm"}, \code{"DV"}, or \code{"SW"}),
#'   calls \code{gradstep()} with the appropriate step-selection method.
#'   Must be length 1 or match \code{length(x)}. Matrices of step sizes are not
#'   supported. Suggestions how to handle all pairs of coordinates are welcome.
#' @param acc.order Integer specifying the desired accuracy order.
#'   The error typically scales as \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.
#' @param side Integer scalar or vector indicating difference type:
#'   \code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
#'   Central differences are recommended unless computational cost is prohibitive.
#' @param symmetric Logical: if \code{TRUE}, then, almost halves computation time
#'   by exploiting Hessian symmetry.
#' @param f0 Optional numeric scalar or vector: if provided and applicable, used
#'   where the stencil contains zero (i.e. \code{FUN(x)} is part of the sum)
#'   to save time.
#'   TODO: Currently ignored.
#' @param h0 Numeric scalar of vector: initial step size for automatic search with
#'   \code{gradstep()}.
#' @param control A named list of tuning parameters passed to \code{gradstep()}.
#' @param cores Integer specifying the number of parallel processes to use. Recommended
#'   value: the number of physical cores on the machine minus one.
#' @param preschedule Logical: if \code{TRUE}, enables pre-scheduling for \code{mclapply()}
#'   and disables load balancing with \code{parLapplyLB()}. Minimises overhead at the cost
#'   of potentially unequal loads at the end of a job. Recommended for functions that
#'   take less than 0.1 s per evaluation.
#' @param func Deprecated; for \code{numDeriv::grad()} compatibility only.
#' @param report Integer: if \code{0}, returns a gradient without any attributes; if \code{1},
#'   attaches the step size and its selection method: \code{2} or higher, attaches the full
#'   diagnostic output (overrides \code{diagnostics = FALSE} in \code{control}).
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @details
#'
#'
#' The optimal step size for 2nd-order-accurate central-differences-based Hessians
#' is of the order Mach.eps^(1/4)
#' to balance the Taylor series truncation error with the rounding error.
#' However, selecting the best step size typically requires knowledge
#' of higher-order cross derivatives and is highly technically involved. Future releases
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
#' @export
#'
#' @seealso [gradstep()] for automatic step-size selection.
#'
#' @examples
#' f <- function(x) prod(sin(x))
#' Hessian(f, 1:4)
#' # Large matrices
#' \donttest{
#'   system.time(Hessian(f, 1:100))
#' }
Hessian <- function(FUN, x, side = 0, acc.order = 2,
                    h = NULL,
                    symmetric = TRUE, h0 = NULL, control = list(), f0 = NULL,
                    cores = 1, preschedule = TRUE, func = NULL, report = 1L, ...) {
  if (is.function(x) && !is.function(FUN)) stop("The argument order must be FUN and then x, not vice versa.")
  n <- length(x)
  h.default <- (abs(x) * (x!=0) + (x==0)) * .Machine$double.eps^(1 / (2 + acc.order))
  if (is.null(h)) h <- h.default
  if (.Platform$OS.type == "windows" && cores > 1) cores <- 1
  ell <- list(...)

  #########################################
  # BEGIN compatibility with numDeriv::grad
  # Detecting numDeriv named arguments (e.g. method.args) in ... first, and handling them
  if (missing(FUN)) {
    if (is.function(func)) {
      FUN <- func
      warning("Use the argument 'FUN' to pass the function for differencing to Hessian instead of 'func'.")
    } else {
      stop("Pass the function for differencing as the named argument 'FUN'.")
    }
  }
  # END compatibility with numDeriv::grad
  #######################################

  if (!is.function(FUN)) stop("'FUN' must be a function.")

  # 'side', 'acc.order', 'h' must align with the length of x
  if (is.null(side)) side <- numeric(n) # NULL --> default central, 0
  if (length(side) == 1) side <- rep(side, n)
  if (!(length(side) %in% c(1, n))) stop("The 'side' argument must have length 1 or same length as x.")
  side[!is.finite(side)] <- 0 # NA --> default 'central -- numDeriv COMPATIBILITY
  if (!all(side %in% -1:1))
    stop("The 'side' argument must contain values 0 for central, 1 for forward, and -1 for backward differences.")

  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)

  # TODO: the part where step is compared to step.CR, step.DV etc.
  if (is.character(h)) {
    stop("Step size algorithms not implemented for Hessians yet.")
  } else if (any(h <= 0)) {
    stop("The argument 'h' (step size) must be positive.")
  }
  if (length(h) == 1) h <- rep(h, n)
  if (length(acc.order) != n) stop("The argument 'acc.order' must have length 1 or length(x).")
  if (length(h) != n) stop("The argument 'h' (step size) must have length 1 or length(x).")

  # TODO: compute f0, check the dimension of f0 output, cannot be a vector
  # TODO: if x is a scalar, do simpler stuff

  xmat.diagonal <- lapply(1:n, function(i) {
    sw <- fdCoef(deriv.order = 2, acc.order = acc.order[i], side = side[i])
    b <- sw$stencil
    w  <- sw$weights
    bh <- sw$stencil * h[i]
    dx <- matrix(0, ncol = n, nrow = length(b))
    dx[, i] <- bh
    xmat <- matrix(rep(x, length(b)), ncol = n, byrow = TRUE) + dx
    if (!is.null(names)) colnames(xmat) <- names(x)
    xmat <- cbind(xmat, index1 = i, index2 = i, weights = w)
    xmat
  })
  xmat.diagonal <- do.call(rbind, xmat.diagonal)
  n1 <- nrow(xmat.diagonal)
  # The duplicated element would be x0
  is.x0 <- which(rownames(xmat.diagonal) == "x")
  is.dup.x0 <- is.x0[-1]

  sw.list <- lapply(1:n, function(j) fdCoef(deriv.order = 1, acc.order = acc.order[j], side = side[j]))
  bh.list <- lapply(1:n, function(j) sw.list[[j]]$stencil * h[j])
  w.list  <- lapply(sw.list, "[[", "weights")
  dx.list <- lapply(1:n, function(j) {
    dx <- matrix(0, ncol = n, nrow = length(bh.list[[j]]))
    dx[, j] <- bh.list[[j]]
    dx
  })
  # TODO: try mixed accuracy orders
  inds <- t(utils::combn(x = n, m = 2))
  colnames(inds) <- c("i", "j")
  # TODO: rewrite this in C++ to eliminate bottlenecks
  xmat.uppertri <- lapply(seq_len(nrow(inds)), function(k) {
    i <- inds[k, "i"]
    j <- inds[k, "j"]
    dxi <- dx.list[[i]]
    dxj <- dx.list[[j]]
    # For each step from dx.list[[i]] and dx.list[[j]]...
    iinds <- cbind(ii = rep(seq_len(nrow(dxi)), nrow(dxj)),
                   jj = rep(seq_len(nrow(dxj)), each = nrow(dxi)))
    xlist <- lapply(seq_len(nrow(iinds)), function(kk) {
      ii <- iinds[kk, "ii"]
      jj <- iinds[kk, "jj"]
      xnew <- x + dxi[ii, ] + dxj[jj, ]
      return(list(x = xnew, weights = w.list[[i]][ii] * w.list[[j]][jj]))
    })
    xmat <- do.call(rbind, lapply(xlist, "[[", "x"))
    xmat <- cbind(xmat, index1 = i, index2 = j,
                  weights = do.call(c, lapply(xlist, "[[", "weights")))
    xmat
  })
  xmat.uppertri <- do.call(rbind, xmat.uppertri)
  n2 <- nrow(xmat.uppertri)

  # De-duplicate from the very beginning, mark f0
  xvals <- rbind(xmat.diagonal, xmat.uppertri)
  rownames(xvals) <- NULL
  weights <- xvals[, "weights"]
  i1 <- as.integer(xvals[, "index1"])
  i2 <- as.integer(xvals[, "index2"])
  spl.ind <- paste0(i1, "_", i2)
  xvals <- xvals[, 1:n, drop = FALSE]
  if (!is.null(names(x))) colnames(xvals) <- names(x)
  xvals <- t(xvals) # Column operations are faster than row operations
  xvals <- lapply(seq_len(ncol(xvals)), function(i) xvals[, i])

  ivals <- c(setdiff(1:n1, is.dup.x0), (n1+1):(n1+n2)) # Unique grid values

  # Parallelising the task in the most efficient way possible, over all values of all grids
  FUN1 <- function(i) do.call(FUN, c(list(x = xvals[[i]]), ell))
  fvals0 <- if (cores > 1) {
    parallel::mclapply(X = ivals, FUN = FUN1, mc.cores = cores, mc.preschedule = preschedule)
  } else {
    lapply(ivals, FUN1)
  }
  if (any(!sapply(fvals0, function(x) is.numeric(x) | is.na(x))))
    stop("'FUN' must output numeric values only, but non-numeric values were returned.")

  # If the output is vector-valued, raise suspicion
  fvals0 <- do.call(rbind, fvals0)
  if (NCOL(fvals0) > 1) stop("Hessian() supports only scalar functions, but the output of FUN has length > 1.")
  fvals0 <- drop(fvals0)
  fvals <- numeric(length(xvals))
  fvals[ivals] <- fvals0
  fvals[is.dup.x0] <- fvals0[is.x0[1]]
  wf <- fvals * weights
  wfh <- wf / unname(h[i1] * h[i2])
  wfh <- split(wfh, f = spl.ind)
  wfh <- vapply(wfh, sum, FUN.VALUE = numeric(1))
  hes <- matrix(NA, n, n)
  hes[lower.tri(hes, diag = TRUE)] <- wfh
  hes[upper.tri(hes)] <- t(hes)[upper.tri(hes)]
  if (!is.null(names(x))) colnames(hes) <- rownames(hes) <- names(x)

  if (report > 0) {
    attr(hes, "step.size") <- h
    # TODO: After implementing autosteps, return the syntax here from Grad
    if (all(h == h.default)) {
      attr(hes, "step.size.method") <- "default"
    } else {
      attr(hes, "step.size.method") <- "user-supplied"
    }
  }

  return(hes)
}
