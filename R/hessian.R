#' Generate grid points for Hessians
#'
#' Creates a list of unique evaluation points for second derivatives: both
#' diagonal (\eqn{\partial^2 / \partial x_i^2}{d^2 / dx^2}) and cross
#' (\eqn{\partial^2 / \partial x_i \partial x_j}{d^2 / dx_i dx_j}).
#'
#' @inheritParams GenD
#'
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{xlist}: a list of unique coordinate shifts,
#'   \item \code{w}: the finite-difference weights (one per point),
#'   \item \code{i1}, \code{i2}: integer vectors giving partial-derivative indices.
#' }
#' The length of each vector matches \code{xlist}.
#' @export
#'
#' @seealso [GenD()], [Hessian()].
#'
#' @examples
#' generateGrid2(1:4, side = rep(0, 4), acc.order = c(2, 6, 4, 2),
#'               h = c(1e-5, 1e-4, 2e-5, 1e-6))
generateGrid2 <- function(x, side, acc.order, h) {
  n <- length(x)

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
  bh.list <- lapply(1:n, function(j) sw.list[[j]]$stencil * h[j])  # Shifts
  w.list  <- lapply(sw.list, "[[", "weights")
  dx.list <- lapply(1:n, function(j) {
    dx <- matrix(0, ncol = n, nrow = length(bh.list[[j]]))
    dx[, j] <- bh.list[[j]]
    dx
  })

  if (n >= 2) {
    inds <- t(utils::combn(x = n, m = 2))
    colnames(inds) <- c("i", "j")

    xmat.offdiag <- lapply(seq_len(nrow(inds)), function(k) {
      i <- inds[k, "i"]
      j <- inds[k, "j"]
      dxi <- dx.list[[i]]
      dxj <- dx.list[[j]]
      # For each step from dx.list[[i]] and dx.list[[j]]...
      iinds <- cbind(ii = rep(seq_len(nrow(dxi)), nrow(dxj)),
                     jj = rep(seq_len(nrow(dxj)), each = nrow(dxi)))
      # Two steps the weights are multiplied
      xlist <- lapply(seq_len(nrow(iinds)), function(kk) {
        ii <- iinds[kk, "ii"]
        jj <- iinds[kk, "jj"]
        return(list(x = x + dxi[ii, ] + dxj[jj, ],
                    weights = w.list[[i]][ii] * w.list[[j]][jj]))
      })
      xmat <- do.call(rbind, lapply(xlist, "[[", "x"))
      xmat <- cbind(xmat, index1 = i, index2 = j,
                    weights = do.call(c, lapply(xlist, "[[", "weights")))
      xmat
    })
    xmat.offdiag <- do.call(rbind, xmat.offdiag)
    n2 <- nrow(xmat.offdiag)
  } else {  # 1x1 Hessians
    xmat.offdiag <- NULL
  }


  # De-duplicate from the very beginning, mark f0
  xvals <- rbind(xmat.diagonal, xmat.offdiag)
  rownames(xvals) <- NULL
  weights <- xvals[, "weights"]
  i1 <- as.integer(xvals[, "index1"])
  i2 <- as.integer(xvals[, "index2"])
  xvals <- xvals[, 1:n, drop = FALSE]  # The true evaluation grid
  if (!is.null(names(x))) colnames(xvals) <- names(x)
  uniq.i <- dupRowInds(xvals)
  xvals <- t(xvals) # Column operations are faster than row operations
  xvals <- lapply(seq_len(ncol(xvals)), function(i) xvals[, i])

  return(list(xlist = xvals, weights = weights, i1 = i1, i2 = i2, uniq.i = uniq.i))
}

#' Numerical Hessians
#'
#' Computes the second derivatives of a function with respect to all combinations
#' of its input coordinates. Arbitrary accuracies and sides for different coordinates
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
#' @param f0 Optional numeric scalar or vector: if provided and applicable, used
#'   where the stencil contains zero (i.e. \code{FUN(x)} is part of the sum)
#'   to save time.
#'   TODO: Currently ignored.
#' @param h0 Numeric scalar of vector: initial step size for automatic search with
#'   \code{gradstep()}.Hessian(f, 1:100)
#' @param control A named list of tuning parameters passed to \code{gradstep()}.
#' @inheritParams runParallel
#' @param func Deprecated; for \code{numDeriv::grad()} compatibility only.
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
#' Some numerical packages use the option (or even the default behaviour) of computing
#' not only the \code{i < j} cross-partials for the Hessian, but all pairs of \code{i}
#' and \code{j}. The upper and lower triangular matrices are filled, and the matrix is
#' averaged with its transpose to obtain a Hessian -- this is the behaviour of \code{optimHess()}.
#' However, it can be shown that \code{H[i, j]} and \code{H[j, i]} use the same evaluation
#' grid, and with a single parallelisable evaluation of the function on that grid, no
#' symmetrisation is necessary because the result is mathematically and computationally identical.
#' In \code{pnd}, only the upper triangular matrix is computed, saving time and ensuring
#' unambiguous results owing to the interchangeability of summation terms (ignoring the numerical
#' error in summation as there is nothing that can be done apart from compensation summation, e.g.
#' via Kahan's algorithm).
#'
#' @return A matrix with as many rows and columns as \code{length(x)}. Unlike the output of
#'   \code{numDeriv::hessian()}, this output preserves the names of \code{x}.
#' @export
#'
#' @seealso [Grad()] for gradients, [GenD()] for generalised numerical differences.
#'
#' @examples
#' f <- function(x) prod(sin(x))
#' Hessian(f, 1:4)
#' # Large matrices
#' \donttest{
#'   system.time(Hessian(f, 1:100))
#' }
Hessian <- function(FUN, x, side = 0, acc.order = 2, h = NULL,
                    h0 = NULL, control = list(), f0 = NULL,
                    cores = 1, preschedule = TRUE, cl = NULL,
                    func = NULL, ...) {
  if (is.function(x) && !is.function(FUN)) {
    warning("The argument order must be FUN and then x, not vice versa.")
    x0 <- FUN
    FUN <- x
    x <- x0
  }

  n <- length(x)

  # 'side', 'acc.order', 'h' must align with the length of x
  if (is.null(side)) side <- numeric(n) # NULL --> default central, 0
  if (length(side) == 1) side <- rep(side, n)
  if (!(length(side) %in% c(1, n))) stop("The 'side' argument must have length 1 or same length as x.")
  side[!is.finite(side)] <- 0 # NA --> default 'central -- numDeriv COMPATIBILITY
  if (!all(side %in% -1:1))
    stop("The 'side' argument must contain values 0 for central, 1 for forward, and -1 for backward differences.")
  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)
  if (length(acc.order) != n) stop("The argument 'acc.order' must have length 1 or length(x).")

  h.default <- stepx(x, deriv.order = 2, acc.order = acc.order)
  if (is.null(h)) h <- h.default

  #########################################
  # BEGIN compatibility with numDeriv::hessian
  ell <- list(...)
  compat <- FALSE
  nd.method <- ell[["method"]]
  nd.method.args <- ell[["method.args"]]
  has.nd.args <- any(names(nd.method.args) %in% c("eps", "d", "zero.tol", "r", "v", "show.details"))
  if ((!is.null(nd.method)) || has.nd.args) {
    compat <- TRUE
    if (is.null(nd.method) && has.nd.args) nd.method <- "Richardson"
    if (length(nd.method) == 1 && nd.method %in% c("simple", "complex", "Richardson")) {
      margs <- ell$method.args
      ma <- list(eps = 1e-4, d = NA, zero.tol = 1e-5, r = 4, show.details = FALSE)
      # Using a slightly smaller step size for one-sided differences
      if (nd.method == "simple") ma$eps <- ma$eps * .Machine$double.eps^(1/4 - 1/5)
      ma[intersect(names(margs), names(ma))] <- margs[intersect(names(margs), names(ma))]
      if (identical(unname(h), unname(h.default))) h <- ma$eps
      if (nd.method == "simple") {
        side <- acc.order <- rep(1L, n)
      } else if (nd.method == "complex") {
        stop("Complex Hessians not implemented yet.")
      } else if (nd.method == "Richardson") {
        side <- numeric(n)
        acc.order <- if (!is.null(ma$r) && is.numeric(ma$r)) 2*ma$r else 8
        if (is.na(ma$d)) ma$d <- .Machine$double.eps^(1 / (2 + acc.order))
        if (is.numeric(ma$v)) {
          warning(paste0("Unlike numDeriv, which uses a large initial step size and ",
                         "shrinkage, pnd uses a smaller initial step and an equispaced ",
                         "symmetric grid. The method argument 'v' will be therefore ignored."))
        }
        is.small <- abs(x) < ma$zero.tol
        h <- (ma$d * abs(x) * !is.small) + (ma$eps * is.small)
      }
      ell[["method"]] <- NULL
    }
    warning(paste0("You are using numDeriv-like syntax. We recommend using the new syntax ",
                   "with more appropriate default values and facilities for automatic ",
                   "step-size selection. See ?Hessian for more information."))
  }

  if (missing(FUN)) {
    if (is.function(func)) {
      FUN <- func
      warning("Use the argument 'FUN' to pass the function for differencing to Hessian instead of 'func'.")
    } else {
      stop("Pass the function for differencing as the named argument 'FUN'.")
    }
  }
  # END compatibility with numDeriv::hessian
  ##########################################

  # Setting up parallel capabilities
  cores <- checkCores(cores)
  if (is.null(cl)) cl <- parallel::getDefaultCluster()

  if (!is.function(FUN)) stop("'FUN' must be a function.")

  # TODO: the part where step is compared to step.CR, step.DV etc.
  if (is.character(h)) {
    stop("Step size algorithms not implemented for Hessians yet.")
  } else if (any(h <= 0)) {
    stop("The argument 'h' (step size) must be positive.")
  }

  if (length(h) == 1) h <- rep(h, n)
  if (length(h) != n) stop("The argument 'h' (step size) must have length 1 or length(x).")

  # TODO: compute f0, check the dimension of f0 output, cannot be a vector
  # TODO: if x is a scalar, do simpler stuff
  # TODO: maybe rewrite this part in C++ to eliminate bottlenecks

  grid <- generateGrid2(x = x, side = side, acc.order = acc.order, h = h)
  ivals <- which(!duplicated(grid$uniq.i))  # Unique grid values
  xvals <- grid$xlist[ivals]  # Unique argument values
  weights <- grid$weights

  # Parallelising the task in the most efficient way possible, over all values of all grids
  fvals <- runParallel(FUN = FUN, x = xvals, cores = cores, cl = cl, preschedule = preschedule)
  fvals <- fvals[grid$uniq.i]  # Restoring duplicates
  nonfinite.f   <- !sapply(fvals, is.finite)
  horrible.f  <- nonfinite.f & (!sapply(fvals, is.na)) & (!sapply(fvals, is.infinite))
  if (any(horrible.f)) {
    warning(paste0("'FUN' must output numeric values only, but some non-numeric values were ",
                   "returned (not even NA or NaN). Some gradient coordinates can be NA. Possible reason: ",
                   "returning character or other type. Check the function output."))
    fvals[horrible.f] <- NA_real_
  } else if (any(nonfinite.f)) {
    warning(paste0("'FUN' must output numeric values only, but some non-numeric values were ",
                   "returned (NA or NaN). Some gradient coordinates can be NA. Possible reason: point at ",
                   "the boundary of the support of FUN. Try side = 1 or -1 for a one-sided solution."))
    fvals[nonfinite.f] <- NA_real_
  }

  # If the output is vector-valued, raise suspicion
  fvals <- do.call(rbind, fvals)
  if (NCOL(fvals) > 1) stop("Hessian() supports only scalar functions, but the output of FUN has length > 1.")
  fvals <- drop(fvals)

  spl.ind <- paste0(grid$i1, "_", grid$i2)
  wf <- fvals * weights
  wfh <- wf / unname(h[grid$i1] * h[grid$i2])
  wfh <- split(wfh, f = spl.ind)
  wfh <- vapply(wfh, sum, FUN.VALUE = numeric(1))
  hes <- matrix(NA, n, n)

  pairs <- strsplit(names(wfh), "_", fixed = TRUE)  # 1_2 becomes c(1, 2)
  ii <- as.integer(vapply(pairs, `[`, 1L, FUN.VALUE = character(1)))
  jj <- as.integer(vapply(pairs, `[`, 2L, FUN.VALUE = character(1)))
  for (k in seq_along(ii)) hes[ii[k], jj[k]] <- wfh[k]
  hes[lower.tri(hes)] <- t(hes)[lower.tri(hes)]

  if (!is.null(names(x))) colnames(hes) <- rownames(hes) <- names(x)

  attr(hes, "step.size") <- h
    # TODO: After implementing autosteps, return the syntax here from Grad
  if (all(h == h.default)) {
    attr(hes, "step.size.method") <- "default"
  } else if (compat) {
    attr(hes, "step.size.method") <- "numDeriv-like"
  } else {
    attr(hes, "step.size.method") <- "user-supplied"
  }

  class(hes) <- "Hessian"
  return(hes)
}
