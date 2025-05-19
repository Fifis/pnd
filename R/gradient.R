#' Determine function dimensionality and vectorisation
#'
#' @inheritParams GenD
#'
#' @details
#' The following combinations of parameters are allowed, suggesting specific input and
#' output handling by other functions:
#'
#' |                               | `elementwise`            | `!elementwise`        |
#' | ----------------------------- | ------------------------ | --------------------- |
#' | `!multivalued`, `vectorised`  | `FUN(xgrid)`             | *(undefined)*         |
#' | `!multivalued`, `!vectorised` | `[mc]lapply(xgrid, FUN)` | `[mc]lapply` gradient |
#' | `multivalued`, `vectorised`   | *(undefined)*            | `FUN(xgrid)` Jacobian |
#' | `multivalued`, `!vectorised`  | *(undefined)*            | `[mc]lapply` Jacobian |
#'
#' Some combinations are impossible: multi-valued functions cannot be element-wise,
#' and single-valued vectorised functions must element-wise.
#'
#' In brief, testing the input and output length and vectorisation capabilities may result in five
#' cases, unlike 3 in \code{numDeriv::grad()} that does not provide checks for Jacobians.
#'
#'
#' @returns A named logical vector indicating if a function is element-wise or not,
#' vectorised or not, and multivalued or not.
#' @order 1
#' @export
#'
#' @examples
#' checkDimensions(sin, x = 1:4, h = 1e-5)  # Rn -> Rn vectorised
#' checkDimensions(function(x) integrate(sin, 0, x)$value, x = 1:4, h = 1e-5)  # non vec
#' checkDimensions(function(x) sum(sin(x)), x = 1:4, h = 1e-5)  # Rn -> R, gradient
#' checkDimensions(function(x) c(sin(x), cos(x)), x = 1, h = 1e-5)  # R -> Rn, Jacobian
#' checkDimensions(function(x) c(sin(x), cos(x)), x = 1:4, h = 1e-5)  # vec Jac
#' checkDimensions(function(x) c(integrate(sin, 0, x)$value, integrate(sin, -x, 0)$value),
#'                  x = 1:4, h = 1e-5)  # non-vectorised Jacobian
checkDimensions <- function(FUN, x, f0 = NULL, func = NULL,
                            elementwise = NA, vectorised = NA, multivalued = NA,
                            deriv.order = 1, acc.order = 2, side = 0, h = NULL,
                            zero.tol = sqrt(.Machine$double.eps), cores = 1, preschedule = TRUE, cl = NULL, ...) {
  if (missing(FUN)) {
    if (is.function(func)) {
      FUN <- func
      warning("Use the argument 'FUN' to pass the function for differencing instead of 'func'.")
    } else {
      stop("Pass the function for differencing as the named argument 'FUN'.")
    }
  }
  if (!is.function(FUN)) stop("'FUN' must be a function.")

  cores <- checkCores(cores)
  if (is.null(cl)) cl <- parallel::getDefaultCluster()

  # Vectorisation checks similar to case1or3 in numDeriv, but better optimised
  # to catch errors and waste fewer evaluations f(x)
  if (!(isTRUE(elementwise) || isFALSE(elementwise) || identical(elementwise, NA)))
    stop("The 'elementwise' argument must be TRUE, FALSE, or NA.")
  if (!(isTRUE(vectorised) || isFALSE(vectorised) || identical(vectorised, NA)))
    stop("The 'vectorised' argument must be TRUE, FALSE, or NA.")
  if (!(isTRUE(multivalued) || isFALSE(multivalued) || identical(multivalued, NA)))
    stop("The 'multivalued' argument must be TRUE, FALSE, or NA.")

  if (isTRUE(multivalued) && isTRUE(elementwise))
    stop("Multi-valued functions cannot be elementwise regardless of vectorisation.")
  if (isFALSE(multivalued) && isFALSE(elementwise) && isTRUE(vectorised))
    stop("Scalar functions with multi-variate inputs cannot be vectorised.")

  # If the user is careful, nothing is necessary
  if ((isTRUE(elementwise) || isFALSE(elementwise)) &&
      (isTRUE(vectorised) || isFALSE(vectorised)) &&
      (isTRUE(multivalued) || isFALSE(multivalued)))
    return(c(elementwise = elementwise, vectorised = vectorised, multivalued = multivalued))

  if (is.null(h) || is.character(h)) {
    h.default <- stepx(x, deriv.order = deriv.order, acc.order = acc.order, zero.tol = zero.tol)
    h <- stats::median(h.default)
  }

  # Making one evaluation and measuring the time
  n <- length(x)
  user.f0 <- !is.null(f0)
  tic0 <- Sys.time()
  if (!user.f0) {  # Try direct evaluation first
    f0 <- safeF(FUN, x, ...)  # If successful, this estimate is almost free from overhead (only CPU spin-up)
  } # Else, we do not know the run time because we are not experimenting with f0
  tic1 <- Sys.time()
  if (checkBadSafeF(f0)) {
    vector.fail <- if (n > 1) TRUE else NA  # Does the function fail on vector input?
    if (is.na(vector.fail))
      stop("Could not evaluate FUN(x) for scalar x. Derivatives or gradients cannot be computed.")
    tic0 <- Sys.time()
    f0 <- suppressWarnings(unlist(runParallel(FUN = function(y) safeF(FUN, y, ...),
                                              x = x, cl = cl, preschedule = preschedule)))
    tic1 <- Sys.time()
  } else {
    vector.fail <- if (n > 1) FALSE else NA
  }
  if (all(!is.finite(f0)))
    warning(paste0("Could not evaluate FUN(x) to at least one finite numeric value ",
                "neither directly nor coordinate-wise. Make sure that the function ",
                "value at the requested point is a numeric vector or scalar."))

  l <- length(f0)

  # The element-wise mapping by definition is having equal lengths of input and output
  # TODO: This is NOT guaranteed, however, to guess correctly:
  # f(x) := c(sum(x), prod(x), sd(x)) evaluated at x = 1:3 would SEEM vectorised; more checks needed

  # 5 valid cases remain after ruling out such things as parallelised + multivalued
  if (n==1 && l==1) {
    elementwise <- TRUE
    multivalued <- FALSE
    vectorised <- NA
  } else if (n>1 && l==1) {  # Gradients
    elementwise <- FALSE
    multivalued <- FALSE
    vectorised <- FALSE
  } else if (n==1 && l>1) {  # Jacobian of a univariate function
    elementwise <- FALSE
    multivalued <- TRUE
    vectorised <- NA
  } else if (n>1 && l>1 && n!=l) {  # Jacobian of a multivariate function
    elementwise <- FALSE
    multivalued <- TRUE
    vectorised <- NA
  } else if (n>1 && l==n) {  # element-wise functions
    elementwise <- TRUE
    multivalued <- FALSE
    vectorised <- !vector.fail
  }

  # TODO: Remaining: what if n = 1, l = 1? Is it vectorised? An extra check is needed!
  # elementwise == TRUE,  multivalued <- FALSE, so we create two extra points: x+h and x-h
  # that can be reused later
  xhvals <- fhvals <- NULL
  if (is.na(vectorised)) { # Works if length(x) == 1, too; creates ((x1-h, ...), (x1+h, ...))
    vectorised <- FALSE

    if (!multivalued) {
      # On a minimal 2-point stencil, check if the output has length 2 and is numeric
      s <- fdCoef(deriv.order = deriv.order[1], acc.order = acc.order[1], side = side[1])
      xhvals  <- x + s$stencil[1:2] * h
      # Guaranteed to safely return NA with an attribute in case of failure
      fhvals <- safeF(FUN, xhvals, ...)
      fherr <- checkBadSafeF(fhvals)
      # If no error and the output has proper length (should be 2), the function handles vector inputs
      if (!fherr && length(fhvals) == length(xhvals)) vectorised <- TRUE
    } else if (n > 1) {  # Can this Jacobian for vectors handle scalars
      fdim <- l / n  # In vectorised Jacobians, the output length should be an integer multiple of n
      xhvals  <- x[1]
      fhvals <- safeF(FUN, xhvals, ...)
      fherr <- checkBadSafeF(fhvals)
      if (!fherr && length(fhvals) == fdim) vectorised <- TRUE  # Integer check is also here
    } else {  # n = 1; can this Jacobian handle vectors
      fdim <- l
      s <- fdCoef(deriv.order = deriv.order[1], acc.order = acc.order[1], side = side[1])
      xhvals  <- x[1] + s$stencil[1:2] * h
      fhvals <- safeF(FUN, xhvals, ...)
      fherr <- checkBadSafeF(fhvals)
      # If no error and the output has proper length (should be 2), the function handles vector inputs
      if (!fherr && length(fhvals) == length(xhvals)*fdim) vectorised <- TRUE
    }
  }

  n0 <- if (n > 1 && n == l) "n" else n
  l0 <- if (l > 1 && n == l) "n" else l
  msg <- paste0("FUN maps R^", n0, " -> R^", l0, " and can", if (!vectorised) "not" else "",
                " process vector inputs. Use 'elementwise = ", elementwise, "', 'vectorised = ",
                vectorised, "', 'multivalued = ", multivalued, "' to skip the checks and save time.")

  ret <- c(elementwise = elementwise, vectorised = vectorised, multivalued = multivalued)
  attr(ret, "x") <- c(x, xhvals)
  if (l>1 && n!=l && length(fhvals)>0) {  # Jacobian checks were made
    fmat <- cbind(matrix(f0, nrow = l, byrow = TRUE), matrix(fhvals, nrow = l, byrow = TRUE))
    attr(ret, "f") <- fmat
  } else {
    attr(ret, "f") <- c(f0, fhvals)
  }
  attr(ret, "seconds") <- if (user.f0) NA else as.numeric(difftime(tic1, tic0, units = "secs"))
  attr(ret, "message") <- msg

  class(ret) <- "checkDimensions"
  return(ret)
}

#' Create a grid of points for a gradient / Jacobian
#'
#' @inheritParams GenD
#' @param stencils A list of outputs from [fdCoef()] for each coordinate of \code{x}.
#'
#' @returns A list with points for evaluation, summation weights for derivative computation, and
#'   indices for combining values.
#' @export
#'
#' @examples
#' generateGrid(1:4, h = 1e-5, elementwise = TRUE, vectorised = TRUE,
#'              stencils = lapply(1:4, function(a) fdCoef(acc.order = a)))
generateGrid <- function(x, h, stencils, elementwise, vectorised) {
  n <- length(x)
  weights <- lapply(stencils, "[[", "weights")
  weights <- unname(unlist(weights))
  stencils <- lapply(stencils, "[[", "stencil")
  slengths <- vapply(stencils, "length", numeric(1))
  index <- rep(1:n, times = slengths)
  if (length(h) == 1) h <- rep(h, n)
  if (length(h) != n) stop("generateGrid: 'h' must have length 1 or length(x).")

  if (elementwise) {  # Apply the stencil to each single coordinate of the input x
    # x is 1-dimensional if the input is repeated scalars
    xlist <- lapply(1:n, function(i)  x[i] + stencils[[i]] * h[i])
    xvals <- do.call(c, xlist)
    if (!is.null(names(x))) names(xvals) <- names(x)[index]
    xvals <- as.list(xvals)
  } else { # Apply the stencil to each coordinate of the full-length x individually
    # for a vector f(x1, x2, ...))
    xlist <- lapply(1:n, function(i) {
      bh <- stencils[[i]] * h[i]
      if (n == 1) {
        xmat <- matrix(x[i] + bh)
      } else {
        dx <- matrix(0, ncol = n, nrow = slengths[i])
        dx[, i] <- bh
        xmat <- matrix(rep(x, slengths[i]), ncol = n, byrow = TRUE) + dx
      }
      xmat
    })
    xvals <- if (vectorised) do.call(c, xlist) else do.call(rbind, xlist)
    if (!is.null(names(x))) colnames(xvals) <- names(x)
    xvals <- t(xvals)  # Column operations are faster than row ones
    xvals <- lapply(seq_len(ncol(xvals)), function(i) xvals[, i])
  }

  return(list(xlist = xvals, weights = weights, index = index))
}

#' Numerical derivative matrices with parallel capabilities
#'
#' Computes numerical derivatives of a scalar or vector function using finite-difference methods.
#' This function serves as a backbone for [Grad()] and [Jacobian()], allowing for detailed control
#' over the derivative computation process, including order of derivatives, accuracy, and step size.
#' \code{GenD} is fully vectorised over different coordinates of the function argument,
#' allowing arbitrary accuracies, sides, and derivative orders for different coordinates.
#'
#' @param FUN A function returning a numeric scalar or a vector whose derivatives are to be
#'   computed. If the function returns a vector, the output will be a Jacobian.
#' @param x Numeric vector or scalar: the point(s) at which the derivative is estimated.
#'   \code{FUN(x)} must be finite.
#' @param elementwise Logical: is the domain effectively 1D, i.e. is this a mapping
#'   \eqn{\mathbb{R} \mapsto \mathbb{R}}{R -> R} or
#'   \eqn{\mathbb{R}^n \mapsto \mathbb{R}^n}{R^n -> R^n}. If \code{NA},
#'   compares the output length ot the input length.
#' @param vectorised Logical: if \code{TRUE}, the function
#'   is assumed to be vectorised: it will accept a vector of parameters and return
#'   a vector of values of the same length. Use \code{FALSE} or \code{"no"}  for
#'   functions that take vector arguments and return outputs of arbitrary length (for
#'   \eqn{\mathbb{R}^n \mapsto \mathbb{R}^k}{R^n -> R^k} functions). If \code{NA},
#'   checks the output length and assumes vectorisation if it matches the input length;
#'   this check is necessary and potentially slow.
#' @param multivalued Logical: if \code{TRUE}, the function is assumed to return vectors longer
#'   than 1. Use \code{FALSE} for element-wise functions. If \code{NA}, attempts inferring it from
#'   the function output.
#' @param deriv.order Integer or vector of integers indicating the desired derivative order,
#'   \eqn{\mathrm{d}^m / \mathrm{d}x^m}{d^m/dx^m}, for each element of \code{x}.
#' @param acc.order Integer or vector of integers specifying the desired accuracy order
#'   for each element of \code{x}.
#'   The final error will be of the order \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.
#' @param stencil Optional custom vector of points for function evaluation.
#'   Must include at least \code{m+1} points for the \code{m}-th order derivative.
#' @param side Integer scalar or vector indicating the type of finite difference:
#'   \code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
#'   Central differences are recommended unless computational cost is prohibitive.
#' @param h Numeric or character specifying the step size(s) for the numerical
#'   difference or a method of automatic step determination (\code{"CR"}, \code{"CRm"},
#'   \code{"DV"}, or \code{"SW"} to be used in [gradstep()]). The default value is
#'   described in \code{?GenD}.
#' @param zero.tol Small positive integer: if \code{abs(x) >= zero.tol}, then, the automatically
#'   guessed step size is relative (\code{x} multiplied by the step), unless an auto-selection
#'   procedure is requested; otherwise, it is absolute.
#' @param f0 Optional numeric: if provided, used to determine the vectorisation type
#'   to save time. If FUN(x) must be evaluated (e.g. second derivatives), saves one evaluation.
#' @param h0 Numeric scalar of vector: initial step size for automatic search with
#'   \code{gradstep()}.
#' @param control A named list of tuning parameters passed to \code{gradstep()}.
#' @inheritParams runParallel
#' @param func For compatibility with \code{numDeriv::grad()} only. If instead of
#'   \code{FUN}, \code{func} is used, it will be reassigned to \code{FUN} with a warning.
#' @param method For compatibility with \code{numDeriv::grad()} only. Supported values:
#'   \code{"simple"} and \code{"Richardson"}. Non-null values result in a warning.
#' @param method.args For compatibility with \code{numDeriv::grad()} only. Check
#'   \code{?numDeriv::grad} for a list of values. Non-empty lists result in a warning.
#' @param ... Additional arguments passed to \code{FUN}. Important! Since R does partial
#'   matching for argument names, using arguments like \code{f(x, a)} and passing \code{a = 1}
#'   to [GenD()], [Grad()], or [Jacobian()] will result in the interpretation \code{acc.order = 1}
#'   because the latter is a named argument of these functions. Either supply the full names of all
#'   similar-looking arguments or name arguments differently (e.g. \code{a0} instead of \code{a}).
#'
#' @details
#'
#'
#' For computing gradients and Jacobians, use convenience wrappers \code{Jacobian} and \code{Grad}.
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
#' The step size  defaults to \code{1.5e-8} (the square root of machine epsilon)
#' for \code{x} less than \code{1.5e-8},
#' \code{(2.2e-16)^(1/(deriv.order + acc.order)) * x} for \code{x > 1}, and interpolates
#' exponentially between these values for \code{1.5e-8 < x < 1}.
#'
#' The use of \code{f0} can reduce computation time similar to the use of \code{f.lower}
#' and \code{f.upper} in \code{uniroot()}.
#'
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
#' @order 1
#' @export
#'
#' @examples
#'
#' # Case 1: Vector argument, vector output
#' f1 <- sin
#' g1 <- GenD(FUN = f1, x = 1:100)
#' g1.true <- cos(1:100)
#' plot(1:100, g1 - g1.true, main = "Approximation error of d/dx sin(x)")
#'
#' # Case 2: Vector argument, scalar result
#' f2 <- function(x) sum(sin(x))
#' g2    <- GenD(FUN = f2, x = 1:4)
#' g2.h2 <- Grad(FUN = f2, x = 1:4, h = 7e-6)
#' g2 - g2.h2  # Tiny differences due to different step sizes
#' g2.auto <- Grad(FUN = f2, x = 1:4, h = "SW")
#' print(attr(g2.auto, "step.search")$exitcode)  # Success
#'
#' # Case 3: vector input, vector argument of different length
#' f3 <- function(x) c(sum(sin(x)), prod(cos(x)))
#' x3 <- 1:3
#' j3 <- GenD(f3, x3, multivalued = TRUE)
#' print(j3)
#'
#' # Gradients for vectorised functions -- e.g. leaky ReLU
#' LReLU <- function(x) ifelse(x > 0, x, 0.01*x)
#' system.time(replicate(10, suppressMessages(GenD(LReLU, runif(30, -1, 1)))))
#' system.time(replicate(10, suppressMessages(GenD(LReLU, runif(30, -1, 1)))))
#'
#' # Saving time for slow functions by using pre-computed values
#' x <- 1:4
#' finner <- function(x) sin(x*log(abs(x)+1))
#' fouter <- function(x) integrate(finner, 0, x, rel.tol = 1e-12, abs.tol = 0)$value
#' # The outer function is non-vectorised
#' fslow <- function(x) {Sys.sleep(0.01); fouter(x)}
#' f0 <- sapply(x, fouter)
#' system.time(GenD(fslow, x, side = 1, acc.order = 2, f0 = f0))
#' # Now, with extra checks, it will be slower
#' system.time(GenD(fslow, x, side = 1, acc.order = 2))
GenD <- function(FUN, x, elementwise = NA, vectorised = NA, multivalued = NA,
                 deriv.order = 1L, side = 0, acc.order = 2L, stencil = NULL,
                 h = NULL, zero.tol = sqrt(.Machine$double.eps),  h0 = NULL, control = list(),
                 f0 = NULL, cores = 1, preschedule = TRUE, cl = NULL,
                 func = NULL, method = NULL, method.args = list(), ...) {
  if (is.function(x) && !is.function(FUN)) {
    warning("The argument order must be FUN and then x, not vice versa.")
    x0 <- FUN
    FUN <- x
    x <- x0
  }

  n <- length(x)

  # 'side', 'deriv.order', 'acc.order', 'h' must align with the length of x
  if (is.null(side)) side <- numeric(n) # NULL --> default central, 0
  if (length(side) == 1) side <- rep(side, n)
  if (!(length(side) %in% c(1, n))) stop("The 'side' argument must have length 1 or length(x).")
  side[!is.finite(side)] <- 0 # NA --> default central for numDeriv compatibility
  if (!all(side %in% -1:1))
    stop("'side' must be 0 for central, 1 for forward, and -1 for backward differences.")
  if (length(deriv.order) == 1) deriv.order <- rep(deriv.order, n)
  if (length(acc.order) == 1) acc.order <- rep(acc.order, n)
  if (length(deriv.order) != n) stop("The argument 'deriv.order' must have length 1 or length(x).")
  if (length(acc.order) != n) stop("The argument 'acc.order' must have length 1 or length(x).")

  h.default <- stepx(x, deriv.order = deriv.order, acc.order = acc.order, zero.tol = zero.tol)
  if (is.null(h)) h <- h.default

  #########################################
  # BEGIN compatibility with numDeriv::grad
  # Detecting numDeriv named arguments in 'method.args'
  compat <- FALSE
  # TODO: check method.args as well
  has.nd.args <- any(names(method.args) %in% c("eps", "d", "zero.tol", "r", "v", "show.details"))
  if ((!is.null(method)) || has.nd.args) {
    compat <- TRUE
    if (is.null(method) && has.nd.args) method <- "Richardson"
    if (length(method) == 1 && method %in% c("simple", "complex", "Richardson")) {
      ma <- list(eps = 1e-5, d = NA, zero.tol = 1e-5, r = 4, show.details = FALSE)
      # Using a better step size for one-sided differences
      if (method == "simple") ma$eps <- ma$eps * .Machine$double.eps^(1/3 - 1/4)
      ma[intersect(names(method.args), names(ma))] <- method.args[intersect(names(method.args), names(ma))]
      if (identical(unname(h), unname(h.default))) h <- ma$eps
      if (method == "simple") {
        side <- acc.order <- rep(1L, n)
      } else if (method == "complex") {
        stop("Complex derivatives not implemented yet.")
      } else if (method == "Richardson") {
        side <- numeric(n)
        acc.order <- if (!is.null(ma$r) && is.numeric(ma$r)) 2*ma$r else 8
        if (is.na(ma$d)) ma$d <- .Machine$double.eps^(1 / (1 + acc.order))
        if (is.numeric(ma$v)) {
          warning(paste0("Unlike numDeriv, which uses a large initial step size and ",
                         "shrinkage, pnd uses a smaller initial step and an equispaced ",
                         "symmetric grid. The method argument 'v' will be therefore ignored."))
        }
        is.small <- abs(x) < ma$zero.tol
        h <- ma$d * abs(x) + ma$eps * is.small
      }
    }
    warning(paste0("You are using numDeriv-like syntax. We recommend using the new syntax ",
                   "with more appropriate default values and facilities for automatic ",
                   "step-size selection. See ?Grad and ?gradstep for more information."))
  }

  if (missing(FUN)) {
    if (is.function(func)) {
      FUN <- func
      warning("Use the argument 'FUN' to pass the function for differencing instead of 'func'.")
    } else {
      stop("Pass the function for differencing as the named argument 'FUN'.")
    }
  }
  # END compatibility with numDeriv::grad
  #######################################

  # Setting up parallel capabilities
  cores <- checkCores(cores)
  if (is.null(cl)) cl <- parallel::getDefaultCluster()

  if (!is.function(FUN)) stop("'FUN' must be a function.")


  # TODO: the part where step is compared to step.CR, step.DV etc.
  # TODO: for long vectorised argument, vectorise the check
  autostep <- FALSE
  if (is.character(h)) {
    hmethod <- h
    h.auto <- gradstep(x = x, FUN = FUN, h0 = h0, method = hmethod, control = control,
                       cores = cores, preschedule = preschedule, cl = cl, ...)
    h <- h.auto$par
    autostep <- TRUE
    # TODO: use this gradient already for central differences
  } else if (any(h <= 0)) {
    stop("The argument 'h' (step size) must be positive.")
  }
  if (length(h) == 1) h <- rep(h, n)
  if (length(h) != n) stop("The argument 'h' (step size) must have length 1 or length(x).")

  # If all deriv.order, acc.order, side are the same, invert Vandermonde matrices only once
  if ((length(unique(deriv.order)) == 1 && length(unique(acc.order)) == 1 &&
      length(unique(side)) == 1) || !is.null(stencil)) {
    # Default stencil or, if not null, the user-supplied one
    stencil1 <- fdCoef(deriv.order = deriv.order[1], acc.order = acc.order[1], side = side[1], stencil = stencil)
    stencils <- replicate(n, stencil1, simplify = FALSE)
  } else {
    stencils <- lapply(1:n, function(i) fdCoef(deriv.order = deriv.order[i],
                                               acc.order = acc.order[i], side = side[i], stencil = stencil))
  }

  # Vectorisation can be inferred if FUN(x) is supplied
  needs.detection <- is.na(elementwise) || is.na(vectorised) || is.na(multivalued)
  if (needs.detection) {
    chk <- checkDimensions(FUN = FUN, x = x, f0 = f0, elementwise = elementwise,
                           vectorised = vectorised, multivalued = multivalued,
                           deriv.order = deriv.order, acc.order = acc.order,
                           side = side, h = h, cl = cl, func = func, zero.tol = zero.tol,
                           cores = cores, preschedule = preschedule, ...)
  } else {
    chk <- c(elementwise = unname(elementwise), vectorised = unname(vectorised),
             multivalued = unname(multivalued))
  }

  grid <- generateGrid(x = x, h = h, stencils = stencils, elementwise = chk["elementwise"],
                       vectorised = chk["vectorised"])

  # Parallelising the task in the most efficient way possible, over all values of all grids
  # TODO: deduplicate, save CPU
  fvals0 <- runParallel(FUN = function(y) safeF(FUN, y, ...), x = grid$x,
                        cores = cores, cl = cl, preschedule = preschedule)
  nonfinite.f   <- !sapply(fvals0, is.finite)
  horrible.f  <- nonfinite.f & (!sapply(fvals0, is.na)) & (!sapply(fvals0, is.infinite))
  if (any(horrible.f)) {
    warning(paste0("'FUN' must output numeric values only, but some non-numeric values were ",
                   "returned (not even NA or NaN). Some gradient coordinates can be NA. Possible reason: ",
                   "returning character or other type. Check the function output."))
    fvals0[horrible.f] <- NA_real_
  } else if (any(nonfinite.f)) {
    warning(paste0("'FUN' must output numeric values only, but some non-numeric values were ",
                   "returned (NA or NaN). Some gradient coordinates can be NA. Possible reason: point at ",
                   "the boundary of the support of FUN. Try side = 1 or -1 for a one-sided solution."))
    fvals0[nonfinite.f] <- NA_real_
  }

  # The function can be vector-valued, which is why we ensure that dimensions are not dropped
  fvals <- do.call(rbind, fvals0)
  wf <- fvals * grid$weights
  wf <- lapply(split(wf, f = grid$index), matrix, ncol = ncol(wf)) # Matrices lose dimensions if split
  wf <- lapply(wf, colSums)
  jac <- unname(do.call(cbind, wf))
  jac <- sweep(jac, 2, h^deriv.order, "/")

  if (!chk["multivalued"]) {
    jac <- drop(jac)
    if (!is.null(names(x))) names(jac) <- names(h) <- names(x)
  } else {
    if (!is.null(names(x))) colnames(jac) <- names(h) <- names(x)
    if (!is.null(colnames(fvals))) rownames(jac) <- colnames(fvals)
  }

  attr(jac, "step.size") <- h
  if (autostep) {
    attr(jac, "step.size.method") <- hmethod
    # TODO: use in printing
    attr(jac, "step.search") <- h.auto
  } else if (all(h == h.default)) {
    attr(jac, "step.size.method") <- "default"
  } else if (compat) {
    attr(jac, "step.size.method") <- "numDeriv-like"
  } else {
    attr(jac, "step.size.method") <- "user-supplied"
  }

  class(jac) <- "GenD"
  return(jac)
}


#' Gradient computation with parallel capabilities
#'
#' Computes numerical derivatives and gradients of scalar-valued functions using
#' finite differences. This function supports both two-sided (central, symmetric) and
#' one-sided (forward or backward) derivatives. It can utilise parallel processing
#' to accelerate computation of gradients for slow functions or
#' to attain higher accuracy faster.
#'
#' @inheritParams GenD
#'
#' @details
#' This function aims to be 100% compatible with the syntax of \code{numDeriv::Grad()},
#' but there might be differences in the step size because some choices made in
#' \code{numDeriv} are not consistent with theory.
#'
#' There is one feature of the default step size in \code{numDeriv} that deserves
#' an explanation. In that package (but not in \code{pnd}),
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
#' @details
#' \code{Grad} does an initial check (if \code{f0 = FUN(x)} is not provided)
#' and calls [GenD()] with a set of appropriate parameters (\code{multivalued = FALSE}
#' if the check succeds). In case of parameter mismatch, throws and error.
#'
#' @return Numeric vector of the gradient. If \code{FUN} returns a vector,
#' a warning is issued suggesting the use of [Jacobian()].
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
#' print(g.auto)
#' attr(g.auto, "step.search")$exitcode  # Success
#'
#' # Gradients for vectorised functions -- e.g. leaky ReLU
#' LReLU <- function(x) ifelse(x > 0, x, 0.01*x)
#' Grad(LReLU, seq(-1, 1, 0.1))
Grad <- function(FUN, x, elementwise = NA, vectorised = NA, multivalued = NA,
                 deriv.order = 1L, side = 0, acc.order = 2, stencil = NULL,
                 h = NULL, zero.tol = sqrt(.Machine$double.eps), h0 = NULL, control = list(),
                 f0 = NULL, cores = 1, preschedule = TRUE, cl = NULL,
                 func = NULL, method = NULL, method.args = list(), ...) {
  if (is.function(x) && !is.function(FUN)) {
    warning("The argument order must be FUN and then x, not vice versa.")
    x0 <- FUN
    FUN <- x
    x <- x0
  }

  cores <- checkCores(cores)
  if (is.null(cl)) cl <- parallel::getDefaultCluster()

  needs.detection <- is.na(elementwise) || is.na(vectorised) || is.na(multivalued)
  if (needs.detection) {
    chk <- checkDimensions(FUN = FUN, x = x, f0 = f0, elementwise = elementwise,
                           vectorised = vectorised, multivalued = multivalued,
                           deriv.order = deriv.order, acc.order = acc.order,
                           side = side, h = h, cl = cl, func = func, zero.tol = zero.tol,
                           cores = cores, preschedule = preschedule, ...)
  } else {
    chk <- c(elementwise = elementwise, vectorised = vectorised, multivalued = multivalued)
  }
  if (chk["multivalued"])
    stop(paste0("Use 'Jacobian()' instead of 'Grad()' for vector-valued functions ",
                "to obtain a matrix of derivatives."))

  d <- GenD(FUN = FUN, x = x,
            elementwise = chk["elementwise"], vectorised = chk["vectorised"], multivalued = chk["multivalued"],
            deriv.order = deriv.order, side = side, acc.order = acc.order, stencil = stencil,
            h = h, zero.tol = zero.tol, h0 = h0, control = control, f0 = f0,
            cores = cores, preschedule = preschedule, cl = cl,
            func = func, method = method, method.args = method.args, ...)

  return(d)
}


#' Jacobian matrix computation with parallel capabilities
#'s
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
#' slowFun <- function(x) {Sys.sleep(0.01); sum(sin(x))}
#' slowFunVec <- function(x) {Sys.sleep(0.01);
#'                            c(sin = sum(sin(x)), exp = sum(exp(x)))}
#' true.g <- cos(1:4)  # Analytical gradient
#' true.j <- rbind(cos(1:4), exp(1:4)) # Analytical Jacobian
#' x0 <- c(each = 1, par = 2, is = 3, named = 4)
#'
#' # Compare computation times
#' system.time(g.slow <- numDeriv::grad(slowFun, x = x0) - true.g)
#' system.time(j.slow <- numDeriv::jacobian(slowFunVec, x = x0) - true.j)
#' system.time(g.fast <- Grad(slowFun, x = x0, cores = 2) - true.g)
#' system.time(j.fast <- Jacobian(slowFunVec, x = x0, cores = 2) - true.j)
#' system.time(j.fast4 <- Jacobian(slowFunVec, x = x0, acc.order = 4, cores = 2) - true.j)
#'
#' # Compare accuracy
#' rownames(j.slow) <- paste0("numDeriv.jac.", c("sin", "exp"))
#' rownames(j.fast) <- paste0("pnd.jac.order2.", rownames(j.fast))
#' rownames(j.fast4) <- paste0("pnd.jac.order4.", rownames(j.fast4))
#' # Discrepancy
#' print(rbind(numDeriv.grad = g.slow, pnd.Grad = g.fast, j.slow, j.fast, j.fast4), 2)
#' # The order-4 derivative is more accurate for functions
#' # with non-zero third and higher derivatives -- look at pnd.jac.order.4
#'
#'
#' @export
Jacobian <- function(FUN, x, elementwise = NA, vectorised = NA, multivalued = NA,
                     deriv.order = 1L, side = 0, acc.order = 2, stencil = NULL,
                     h = NULL, zero.tol = sqrt(.Machine$double.eps), h0 = NULL,
                     control = list(), f0 = NULL,
                     cores = 1, preschedule = TRUE, cl = NULL,
                     func = NULL, method = NULL, method.args = list(), ...) {
  if (is.function(x) && !is.function(FUN)) {
    warning("The argument order must be FUN and then x, not vice versa.")
    x0 <- FUN
    FUN <- x
    x <- x0
  }

  cores <- checkCores(cores)
  if (is.null(cl)) cl <- parallel::getDefaultCluster()

  needs.detection <- is.na(elementwise) || is.na(vectorised) || is.na(multivalued)
  if (needs.detection) {
    chk <- checkDimensions(FUN = FUN, x = x, f0 = f0,
                           elementwise = elementwise, vectorised = vectorised, multivalued = multivalued,
                           deriv.order = deriv.order, acc.order = acc.order,
                           side = side, h = h, cl = cl, func = func, zero.tol = zero.tol,
                           cores = cores, preschedule = preschedule, ...)
  } else {
    chk <- c(elementwise = elementwise, vectorised = vectorised, multivalued = multivalued)
  }
  # if (!chk["multivalued"])
  #   warning(paste0("Use 'Grad()' instead of 'Jacobian()' for scalar-valued functions ",
  #               "to obtain a vector of derivatives."))
  # Instead of this old check honour the user's request and really return
  # a row vector for 1D output

  d <- GenD(FUN = FUN, x = x,
            elementwise = chk["elementwise"], vectorised = chk["vectorised"], multivalued = chk["multivalued"],
            deriv.order = deriv.order, side = side, acc.order = acc.order, stencil = stencil,
            h = h, h0 = h0, zero.tol = zero.tol, control = control, f0 = f0,
            cores = cores, preschedule = preschedule, cl = cl,
            func = func, method = method, method.args = method.args, ...)

  if (is.null(dim(d))) {
    a <- attributes(d)
    l <- length(d)
    nms <- a[["names"]]
    a$names <- NULL
    d <- matrix(d, nrow = 1)
    attributes(d)[names(a)] <- a
    colnames(d) <- nms
  }

  return(d)
}
