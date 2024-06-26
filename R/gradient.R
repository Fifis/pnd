#' Numerical derivative matrices with parallel capabilities
#'
#' Computes numerical derivatives of a scalar or vector function using finite-difference methods.
#' This function serves as a backbone for [Grad()] and [Jacobian()], allowing for detailed control
#' over the derivative computation process, including order of derivatives, accuracy, and step size.
#' \code{GenD} is fully vectorised over different coordinates of the function argument,
#' allowing arbitrary accuracies, sides, and derivative orders for different coordinates.
#'
#' @param FUN A function returning a numeric scalar or a vector whose derivatives are to be computed
#'   If the function returns a vector, the output will be is a Jacobian.
#' @param x Numeric vector or scalar: the point(s) at which the derivative is estimated.
#'   \code{FUN(x)} must be finite value.
#' @param deriv.order Integer or vector of integers indicating the desired derivative order,
#'   \eqn{\mathrm{d}^m / \mathrm{d}x^m}{d^m/dx^m}, for each element of \code{x}.
#' @param acc.order Integer or vector of integers specifying the desired accuracy order
#'   for each element of \code{x}.
#'   The final error will be of the order \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.
#' @param side Integer scalar or vector indicating the type of finite difference:
#'   \code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
#'   Central differences are recommended unless computational cost is prohibitive.
#' @param h Numeric or character specifying the step size(s) for the numerical
#'   difference or a method of automatic step determination (\code{"CR"}, \code{"CRm"},
#'   \code{"DV"}, or \code{"SW"} to be used in [gradstep()]).
#' @param f0 Optional numeric: if provided and applicable, used
#'   where the stencil contains zero (i.e. \code{FUN(x)} is part of the sum)
#'   to save time.
#'   TODO: Currently ignored.
#' @param h0 Numeric scalar of vector: initial step size for automatic search with
#'   \code{gradstep()}.
#' @param control A named list of tuning parameters passed to \code{gradstep()}.
#' @param cores Integer specifying the number of CPU cores used for parallel computation.
#'   Recommended to be set to the number of physical cores on the machine minus one.
#' @param load.balance Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
#'   or enables load balancing with \code{parLapplyLB()}.
#' @param func For compatibility with \code{numDeriv::grad()} only. If instead of
#'   \code{FUN}, \code{func} is used, it will be reassigned to \code{FUN} with a warning.
#' @param report Integer for the level of detail in the output. If \code{0},
#'   returns a gradient without any attributes; if \code{1},
#'   attaches the step size and its selection method: \code{2} or higher attaches the full
#'   diagnostic output as an attribute.
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @details
#'
#'
#' For computation of Jacobians, use \code{Jacobian} or \code{Grad}. These two functions
#' are equivalent, but using \code{Grad} for vector-valued returns will produce a warning.
#'
#' If the step size is too large, the slope of the secant poorly estimates the derivative;
#' if it is too small, it leads to numerical instability due to the function value rounding.
#'
#' The optimal step size for one-sided differences typically approaches Mach.eps^(1/2)
#' to balance the Taylor series truncation error with the rounding error due to storing
#' function values with limited precision. For two-sided differences, it is proportional
#' to Mach.eps^(1/3). However, selecting the best step size typically requires knowledge
#' of higher-order derivatives, which may not be readily available. Luckily,
#' using \code{step = "SW"} invokes a reliable automatic data-driven step-size selection.
#' Other options include \code{"DV"}, \code{"CR"}, and \code{"CRm"}.
#'
#' The use of \code{f0} can reduce computation time similar to the use of \code{f.lower}
#' and \code{f.upper} in \code{uniroot()}.
#'
#' For convenience, \code{report = 2} overrides \code{diagnostics = FALSE} in the
#' \code{control}) list.
#'
#' Unlike \code{numDeriv::grad()} and \code{numDeriv::jacobian()}, this function
#' fully preserves the names of \code{x} and \code{FUN(x)}.
#'
#' @return A vector or matrix containing the computed derivatives, structured according
#'   to the dimensionality of \code{x} and \code{FUN}. If \code{FUN} is scalar-valued,
#'   returns a gradient vector. If \code{FUN} is vector-valued, returns a Jacobian matrix.
#'
#' @seealso [gradstep()] for automatic step-size selection.
#'
#' @export
GenD <- function(FUN, x,
                 deriv.order = 1L, side = 0, acc.order = 2L,
                 h = (abs(x)*(x!=0) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order)),
                 h0 = NULL, control = list(), f0 = NULL,
                 cores = 1, load.balance = TRUE, func = NULL,
                 report = 1L,
                 ...) {
  if (is.function(x) && !is.function(FUN)) stop("The argument order must be FUN and then x, not vice versa.")
  n <- length(x)
  if (.Platform$OS.type == "windows" && cores > 1) cores <- 1
  h.default <- (abs(x) * (x!=0) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order))
  #########################################
  # BEGIN compatibility with numDeriv::grad
  # Detecting numDeriv named arguments (e.g. method.args) in ... first, and handling them
  ell <- list(...)
  compat <- FALSE
  if (!is.null(m <- ell[["method"]])) {
    compat <- TRUE
    if (length(m) == 1 && m %in% c("simple", "complex", "Richardson")) {
      margs <- ell$method.arguments
      ma <- list(eps = 1e-5, d = NA, zero.tol = 1e-5, r = 4, show.details = FALSE)
      # Using a better step size for one-sided differences
      if (m == "simple") ma$eps <- (abs(x) * (x!=0) + (x==0)) * sqrt(.Machine$double.eps) * 2
      ma[names(margs)] <- margs
      # If the user did not supply a custom step size, use numDeriv method arguments
      if (identical(unname(h), unname(h.default))) h <- ma$eps
      if (m == "simple") {
        side <- acc.order <- rep(1L, n)
      } else if (m == "complex") {
        stop("Complex derivatives not implemented yet.")
      } else if (m == "Richardson") {
        side <- numeric(n)
        acc.order <- if (!is.null(ma$r) && is.numeric(ma$r)) 2*ma$r else 8
        if (is.na(ma$d)) ma$d <- .Machine$double.eps^(1 / (1 + acc.order))
        if (is.numeric(ma$v)) {
          warning(paste0("Unlike numDeriv, which uses a large initial step size and ",
                         "shrinkage, pnd uses a smaller initial step and an enlarging equispaced ",
                         "grid. The method argument 'v' will be therefore ignored."))
        }
        is.small <- abs(x) < ma$zero.tol
        h <- ma$d * abs(x) + ma$eps * is.small
      }
      ell[["method"]] <- NULL
    }
    warning(paste0("You are using numDeriv-like syntax. We recommend using the new syntax ",
                   "with more approproate default values and facilities for automatic ",
                   "step-size selection. See ?Grad and ?gradstep for more information."))
  }

  if (missing(FUN)) {
    if (is.function(func)) {
      FUN <- func
      warning("Use the argument 'FUN' to pass the function for differencing to Grad instead of 'func'.")
    } else {
      stop("Pass the function for differencing as the named argument 'FUN'.")
    }
  }
  # END compatibility with numDeriv::grad
  #######################################

  if (!is.function(FUN)) stop("'FUN' must be a function.")

  # 'side', 'deriv.order', 'acc.order', 'h' must align with the length of x
  if (is.null(side)) side <- numeric(n) # NULL --> default central, 0
  if (length(side) == 1) side <- rep(side, n)
  if (!(length(side) %in% c(1, n))) stop("The 'side' argument must have length 1 or same length as x.")
  side[!is.finite(side)] <- 0 # NA --> default 'central -- numDeriv COMPATIBILITY
  if (!all(side %in% -1:1)) stop("The 'side' argument must contain values 0 for central, 1 for forward, and -1 for backward difference.")

  if (length(deriv.order) == 1) deriv.order <- rep(deriv.order, n)
  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)

  # TODO: the part where step is compared to step.CR, step.DV etc.
  autostep <- FALSE
  if (is.character(h)) {
    method <- h
    if (report == 2) control$diagnostics <- TRUE
    h.auto <- gradstep(x = x, FUN = FUN, h0 = h0, method = method, control = control, ...)
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
  # The rubbish in the ellipsis (...) should have been sanitised by now
  FUN1 <- function(x) do.call(FUN, c(list(x = x), ell))
  fvals <- if (cores > 1) {
    parallel::mclapply(X = xvals, FUN = FUN1, mc.cores = cores, mc.preschedule = !load.balance)
  } else {
    lapply(xvals, FUN1)
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
    if (!is.null(names(x))) names(jac) <- names(h) <- names(x)
  } else {
    if (!is.null(names(x))) colnames(jac) <- names(h) <- names(x)
    if (!is.null(colnames(fvals))) rownames(jac) <- colnames(fvals)
  }

  if (report > 0) {
    attr(jac, "step.size") <- h
    if (autostep) {
      attr(jac, "step.size.method") <- method
    } else if (all(h == (abs(x) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order)))) {
      attr(jac, "step.size.method") <- "default"
    } else if (compat) {
      attr(jac, "step.size.method") <- "numDeriv-like"
    } else {
      attr(jac, "step.size.method") <- "user-supplied"
    }
    if (autostep && report > 1) attr(jac, "step.search") <- h.auto
  }

  return(jac)
}


#' Gradient computation with parallel capabilities
#'
#' Computes numerical derivatives and gradients of scalar-valued functions using
#' finite differences. This function supports both two-sided (central, symmetric) and
#' one-sided (forward or backward) derivatives. It can utilise parallel processing
#' to accelerate computation of gradients for slow functions or
#' to attain higher accuracy faster. Currently, only Mac and Linux are supported
#' \code{parallel::mclapply()}. Windows support with \code{parallel::parLapply()}
#' is under development.
#'
#' @inheritParams GenD
#'
#' @details
#' This function aims to be 100% compatible with the syntax of \code{numDeriv::Grad()}.
#'
#' There is one feature of the default step size in \code{numDeriv} that deserves
#' an explanation.
#'
#' \itemize{
#'   \item If \code{method = "simple"}, then, simple forward differences are used with
#'   a fixed step size \code{eps}, which we denote by \eqn{\varepsilon}{eps}.
#'   \item If \code{method = "Richardson"}, then, central differences are used with
#'   a fixed step
#'   \eqn{h := |d\cdot x| + \varepsilon (|x| < \mathrm{zero.tol})}{h := |d*x| + eps*(|x| < zero.tol)},
#'   where \code{d = 1e-4} is the relative step size and \code{eps} becomes an extra
#'   addition to the step size for the argument that are closer to zero than \code{zero.tol}.
#' }
#' We believe that the latter may lead to mistakes when the user believes that they can set
#' the step size for near-zero arguments, whereas in reality, a combination of \code{d} and \code{eps}
#' is used.
#'
#' Here is the synopsis of the old arguments:
#' \describe{
#'   \item{side}{\code{numDeriv} uses \code{NA} for handling two-sided differences.
#'   The \code{pnd} equivalent is \code{0}, and \code{NA} is replaced with \code{0}.}
#'   \item{eps}{If \code{numDeriv} \code{method = "simple"}, then, \code{eps = 1e-4} is
#'   the absolute step size and forward differences are used.
#'   If \code{method = "Richardson"}, then, \code{eps = 1e-4} is the absolute increment of the step
#'   size for small arguments below the zero tolerance.}
#'   \item{d}{If \code{numDeriv} \code{method = "Richardson"}, then, \code{d*abs(x)} is the
#'   step size for arguments above the zero tolerance and the baseline step size for
#'   small arguments that gets incremented by \code{eps}.}
#'   \item{r}{The number of Richardson extrapolations that successively reduce the initial step size.
#'   For two-sided differences, each extrapolation increases the accuracy order by 2.}
#'   \item{v}{The reduction factor in Richardson extrapolations.}
#' }
#'
#' Here are the differences in the new compatible implementation.
#' \describe{
#'   \item{eps}{If \code{numDeriv} \code{method = "simple"}, then,
#'   \code{ifelse(x!=0, abs(x), 1) * sqrt(.Machine$double.eps) * 2} is used because
#'   one-sided differences require a smaller step size to reduce the truncation error.
#'   If \code{method = "Richardson"}, then, \code{eps = 1e-5}.}
#'   \item{d}{If \code{numDeriv} \code{method = "Richardson"}, then, \code{d*abs(x)} is the
#'   step size for arguments above the zero tolerance and the baseline step size for
#'   small arguments that gets incremented by \code{eps}.}
#'   \item{r}{The number of Richardson extrapolations that successively reduce the initial step size.
#'   For two-sided differences, each extrapolation increases the accuracy order by 2.}
#'   \item{v}{The reduction factor in Richardson extrapolations.}
#' }
#'
#' @return Numeric vector of the gradient. If \code{FUN} returns a vector,
#' a warning is issued suggesting the use of `Jacobian()`.
#'
#' @seealso [GenD()], [Jacobian()]
#' @export
#'
#' @examples
#' f <- function(x) sum(sin(x))
#' g1 <- Grad(FUN = f, x = 1:4)
#' g2 <- Grad(FUN = f, x = 1:4, h = 7e-6)
#' g2 - g1  # Tiny differences due to different step sizes
#' g.auto <- Grad(FUN = f, x = 1:4, h = "SW")
#' g3.full <- Grad(FUN = f, x = 1:4, h = "SW", report = 2)
#' print(g3.full)
#' attr(g3.full, "step.search")$exitcode  # Success
#'
Grad <- function(FUN, x, deriv.order = 1L, side = 0, acc.order = 2,
                 h = (abs(x)*(x!=0) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order)),
                 h0 = NULL, control = list(), f0 = NULL,
                 cores = 1, load.balance = TRUE,
                 func = NULL, report = 1L, ...) {
  d <- GenD(FUN = FUN, x = x, deriv.order = deriv.order, side = side, acc.order = acc.order,
            h = h, h0 = h0, control = control, f0 = f0, cores = cores, load.balance = load.balance,
            func = func, report = report, ...)
  if (is.matrix(d))
    warning(paste0("Use 'Jacobian()' instead of 'Grad()' for vector-valued functions ",
                   "to obtain a matrix of derivatives."))
  return(d)
}


#' Jacobian matrix computation with parallel capabilities
#'
#' Computes the numerical Jacobian for vector-valued functions. Its columns are
#' partial derivatives of the function with respect to the input elements.
#' This function supports both two-sided (central, symmetric) and
#' one-sided (forward or backward) derivatives. It can utilise parallel processing
#' to accelerate computation of gradients for slow functions or
#' to attain higher accuracy faster. Currently, only Mac and Linux are supported
#' \code{parallel::mclapply()}. Windows support with \code{parallel::parLapply()}
#' is under development.
#'
#' @inheritParams GenD
#'
#' @return Matrix where each row corresponds to a function output and each column
#' to an input coordinate. For scalar-valued functions, a warning is issued and
#' the output is returned as a row matrix.
#'
#' @seealso [GenD()], [Grad()]
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
Jacobian <- function(FUN, x,
                     deriv.order = 1L, side = 0, acc.order = 2,
                     h = (abs(x)*(x!=0) + (x==0)) * .Machine$double.eps^(1 / (deriv.order + acc.order)),
                     h0 = NULL, control = list(), f0 = NULL,
                     cores = 1, load.balance = TRUE, func = NULL,
                     report = 1L, ...) {
  d <- GenD(FUN = FUN, x = x, deriv.order = deriv.order, side = side, acc.order = acc.order,
            h = h, h0 = h0, control = control, f0 = f0, cores = cores, load.balance = load.balance,
            func = func, report = report, ...)
  if (is.null(dim(d))) {
    warning(paste0("Use 'Grad()' instead of 'Jacobian()' for scalar-valued functions ",
                   "to obtain a vector of derivatives. This output is a matrix with 1 row, ",
                   "but a vector with NULL dimensions would be more appropriate."))
    d <- matrix(d, nrow = 1)
  }
  return(d)
}
