#' Numerically stable non-confluent Vandermonde system solver
#'
#' @param s Numeric vector of stencil points defining the Vandermonde matrix on
#'   the left-hand side, where each element \eqn{S_{i, j}}{S_ij} is calculated as
#'   \code{s[j]^(i-1)}.
#' @param b Numeric vector of the right-hand side of the equation.
#'   This vector must be the same length as \code{s}.
#'
#' @details
#' This function utilises the \insertCite{bjorck1970solution}{pnd} algorithm for
#' an accurate solution to non-confluent Vandermonde systems,
#' which are known for their numerical instability. Unlike Gaussian elimination,
#' which suffers from ill conditioning, this algorithm achieves numerical stability
#' through exploiting the ordering of the stencil.
#' An unsorted stencils will trigger a warning.
#' Additionally, the stencil must contain unique points, as repeated values make
#' the Vandermonde matrix confluent and therefore non-invertible.
#'
#' This implementation is a verbatim translation of Algorithm 4.6.2 from
#' \insertCite{golub2013matrix}{pnd}, which is robust against the issues
#' typically associated with Vandermonde systems.
#'
#' See \insertCite{higham1987error}{pnd} for an in-depth error analysis of this algorithm.
#'
#' @return A numeric vector of coefficients solving the Vandermonde system,
#'   matching the length of \code{s}.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' # Approximate the 4th derivatives on a non-negative stencil
#' solveVandermonde(s = 0:5, b = c(0, 0, 0, 0, 24, 0))
#'
#' # Small numerical inaccuracies: note the 6.66e-15 in the 4th position --
#' # it should be rounded towards zero:
#' solveVandermonde(s = -3:3, b = c(0, 1, rep(0, 5))) * 60
solveVandermonde <- function(s, b) {
  if (length(unique(s)) != length(s)) stop("The stencil points in 's' must be unique.")
  if (length(s) != length(b)) stop("The input arguments 's' and 'b' must be of the same length.")
  if (is.unsorted(s)) warning(paste0("Numerical solutions of Vandermonde systems based on ",
    "unsorted inputs are EXTREMELY unstable. Sort the input in ascending order!"))
  # Algorithm 4.6.2 from Golub & Loan (2013) 'Matrix computations', 4th ed.
  n <- length(s) - 1
  for (k in 0:(n-1)) {
    for (i in n:(k+1)) b[i+1] <- b[i+1] - s[k+1] * b[i]
  }
  for (k in (n-1):0) {
    for (i in (k+1):n) b[i+1] <- b[i+1] / (s[i+1] - s[i-k])
    for (i in k:(n-1)) b[i+1] <- b[i+1] - b[i+2]
  }
  return(b)
}


#' Finite-difference coefficients for arbitrary grids
#'
#' This function computes the coefficients for a numerical approximation to derivatives
#' of any specified order. It provides the minimally sufficient stencil
#' for the chosen derivative order and desired accuracy order.
#' It can also use any user-supplied stencil (uniform or non-uniform).
#' For that stencil \eqn{\{b_i\}_{i=1}^n}{{b[i]}, i = 1, ..., n}, it computes
#' the optimal weights \eqn{\{w_i\}}{{w[i]}} that yield
#' the numerical approximation of the derivative:
#' \deqn{\frac{d^m f}{dx^m} \approx h^{-m} \sum_{i=1}^n w_i f(x + b_i\cdot h)}{d^m/dx^m f(x) ~ sum_i w[i] f(x + b[i]*h)}
#'
#' @param deriv.order Order of the derivative (\eqn{m}{m} in \eqn{\frac{d^m f}{dx^m}}{d^m/dx^m f(x)})
#'   for which a numerical approximation is needed.
#' @param acc.order Order of accuracy: defines how the approximation error scales
#'   with the step size \eqn{h}{h}, specifically \eqn{O(h^{a+1})}{O(h^(a+1))}, where
#'   \eqn{a}{a} is the accuracy order and depends on the higher-order derivatives
#'   of the function.
#' @param side Integer that determines the type of finite-difference scheme:
#'   \code{0} for central (AKA symmetrical or two-sided; the default),
#'   \code{1} for forward, and \code{-1} for backward.
#'   Using \code{2} (for 'two-sided') triggers a warning and is treated as \code{0}.
#'   with a warning. Unless the function is computationally prohibitively,
#'   central differences are strongly recommended for their accuracy.
#' @param stencil Optional custom vector of points for function evaluation.
#'   Must include at least \code{m+1} points for the \code{m}-th order derivative.
#' @param zero.action Character string specifying how to handle near-zero weights:
#'   \code{"drop"} to omit small (less in absolute value than \code{zero.tol} times
#'   the median weight) weights and corresponding stencil points, \code{"round"}
#'   to round small weights to zero, and \code{"none"} to leave
#'   all weights as calculated.
#'   E.g. the stencil for \eqn{f'(x)}{f'(x)} is \code{(-1, 0, 1)} with weights
#'   \code{(-0.5, 0, 0.5)}; using \code{"drop"} eliminates the zero weight,
#'   and the redundant \code{f(x)} is not computed.
#' @param zero.tol Non-negative scalar defining the threshold: weights below
#'   \code{zero.tol} times the median weight are considered near-zero.
#'
#' @details
#' This function relies on the approach of approximating numerical derivarives
#' by weghted sums of function values described in \insertCite{fornberg1988generation}{pnd}.
#' It reproduces all tables from this paper exactly; see the example below to
#' create Table 1.
#'
#' The finite-difference coefficients for any given stencil are given as a solution of a linear
#' system. The capabilities of this function are similar to those of \insertCite{taylor2016finite}{pnd},
#' but instead of matrix inversion, the \insertCite{bjorck1970solution}{pnd} algorithm is used because
#' the left-hand-side matrix is a Vandermonde matrix, and its inverse may be
#' very inaccurate, especially for long one-sided stencils.
#'
#' The weights computed for the stencil via this algorithm are very reliable; numerical
#' simulations in \insertCite{higham1987error}{pnd} show that the relative error is
#' low even for ill-conditioned systems. \insertCite{kostyrka2024what}{pnd}
#' computes the exact relative error of the weights on the stencils returned by
#' this function; the zero tolerance is based on these calculations.
#'
#' @return A list containing the \code{stencil} used and the corresponding
#'   \code{weights} for each point.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' fdCoef()  # Simple two-sided derivative
#' fdCoef(2) # Simple two-sided second derivative
#' fdCoef(acc.order = 4)$weights * 12  # Should be (1, -8, 8, -1)
#'
#' # Using an custom stencil for the first derivative: x-2h and x+h
#' fdCoef(stencil = c(-2, 1), acc.order = 1)
#'
#' # Reproducing Table 1 from Fornberg (1988) (cited above)
#' pad9 <- function(x) {l <- length(x); c(a <- rep(0, (9-l)/2), x, a)}
#' f <- function(d, a) pad9(fdCoef(deriv.order = d, acc.order = a,
#'                                 zero.action = "round")$weights)
#' t11 <- t(sapply((1:4)*2, function(a) f(d = 1, a)))
#' t12 <- t(sapply((1:4)*2, function(a) f(d = 2, a)))
#' t13 <- t(sapply((1:3)*2, function(a) f(d = 3, a)))
#' t14 <- t(sapply((1:3)*2, function(a) f(d = 4, a)))
#' t11 <- cbind(t11[, 1:4], 0, t11[, 5:8])
#' t13 <- cbind(t13[, 1:4], 0, t13[, 5:8])
#' t1 <- data.frame(OrdDer = rep(1:4, times = c(4, 4, 3, 3)),
#'                  OrdAcc = c((1:4)*2, (1:4)*2, (1:3)*2, (1:3)*2),
#'                  rbind(t11, t12, t13, t14))
#' colnames(t1)[3:11] <- as.character(-4:4)
#' print(t1, digits = 4)
fdCoef <- function(deriv.order = 1L, side = c(0L, 1L, -1L),
                   acc.order = 2L, stencil = NULL,
                   zero.action = c("drop", "round", "none"),
                   zero.tol = NULL) {
  if (length(deriv.order) != 1L) stop("The 'deriv.order' argument must be an integer of length 1.")
  if (length(acc.order) != 1L) stop("The 'acc.order' argument must must be an integer of length 1.")
  if (!is.numeric(deriv.order) || deriv.order != round(deriv.order) || deriv.order < 0)
    stop("'deriv.order' must be a non-negative integer of length 1.")
  if (!is.numeric(acc.order) || acc.order != round(acc.order) || acc.order <= 0)
    stop("'acc.order' must be a positive integer of length 1.")
  side <- side[1]
  if (side == 2) {
    side <- 0 # Interpreting '2' as two-sided = central
    warning("Interpreting 'side = 2' as 'two-sided central differences'; please use side = 0.")
  }
  if (!(side %in% -1:1))
    stop("The 'side' argument must be -1, 0, or 1 for backward, central, and forward differences.")
  zero.action <- zero.action[1]
  if (!(zero.action %in% c("drop", "round", "none")))
    stop("The 'zero.action' argument must be 'drop', 'round', or 'none'.")
  # TODO: implement interpolation
  if (deriv.order == 0) stop("Derivative order 0 (interpolation) not implemented yet.")

  if (is.null(stencil)) { # Miminally sufficient stencil to prevent zero weights
    acc.requested <- acc.order
    if (side == -1L) {    # and exponential growth of the Vandermonde matrix elements
      stencil <- (-acc.order - deriv.order + 1):0L
    } else if (side == 1L) {
      stencil <- 0L:(acc.order + deriv.order - 1)
    } else { # Central differences
      if (acc.order %% 2 == 1) {
        prefix <- switch(acc.order, "st", "nd", "rd")
        if (is.null(prefix)) prefix <- "th"
        acc.order <- acc.order + 1L
        warning(paste0("You requested ", acc.order-1, prefix, "-order-accurate central differences, ",
                       "but the minimal stencil will provide accuracy order ", acc.order,
                       " -- higher by 1, which is generally desirable."))
      }
      l.end <- floor(acc.order/2) + floor((deriv.order-1)/2)
      stencil <- (-l.end):l.end
      if (deriv.order %% 2 == 1) stencil <- setdiff(stencil, 0)
    }
  } else {
    acc.requested <- NA
  }

  stencil <- sort(stencil)
  is.dup <- duplicated(stencil)
  if (any(is.dup)) {
    warning(paste0("The user-supplied stencil contains duplicates: ",
                   paste0(stencil[is.dup], collapse = ", "), " -- dropping them."))
    stencil <- stencil[!is.dup]
  }
  l <- length(stencil)
  is.symm <- all(stencil + rev(stencil) == 0)
  a <- l - deriv.order + is.symm # Effective accuracy (preliminary guess because 0 may be redundant)
  if (l < deriv.order + 1)
    stop("To compute the m-th derivative, at least m+1 unique stencil points are required.")
  if (l < deriv.order + acc.order - is.symm) {
    prefix <- switch(a, "st", "nd", "rd")
    if (is.null(prefix)) prefix <- "th"
    warning(paste0("The user-supplied stencil needs ", acc.order - a,
                   " more unique points to achieve the requested accuracy order ", acc.order, ". ",
                   "The result will be only ", a, prefix, "-order-accurate."))
  }

  b <- numeric(l)
  b[1 + deriv.order] <- factorial(deriv.order)
  weights <- solveVandermonde(s = stencil, b = b)

  if (is.null(zero.tol)) {
    rel.tol <- 1024 * .Machine$double.eps
    if (a > 6) rel.tol <- rel.tol * 10^((min(a, 16)-6)/2)
    zero.tol <- max(stats::median(abs(weights)), 1) * rel.tol
  }
  if (zero.action != "none") {
    zw <- abs(weights) < zero.tol
    if (zero.action == "drop") {
      stencil <- stencil[!zw]
      weights <- weights[!zw]
    } else if (zero.action == "round") {
      weights[zw] <- 0
    }
  }
  rs <- round(stencil, 2) # Suitable names for reasonable integer stencils
  if (any(duplicated(rs))) {
    rs <- sprintf("%1.2e", stencil)
    names(weights) <- paste0("x", ifelse(stencil < 0, "", "+"), rs, "h")
    warning(paste("Stencils should contain large entries like '-0.5, 1, 3', but this one has",
                  "two very close points. Wildly uneven gaps may result in numerical instability."))
  } else {
    names(weights) <- paste0("x", ifelse(stencil < 0, "-", "+"), abs(rs), "h")
  }
  names(weights)[stencil == 0] <- "x"

  # Computing the coefficient on the remainder term in the output
  oseq <- 0:(deriv.order+a)
  B <- outer(stencil, oseq, "^")
  resulting.terms <- colSums(B * weights) / factorial(oseq)
  resulting.terms[abs(resulting.terms) < zero.tol] <- 0
  flabs <- sapply(oseq, function(i) if (i == 0) " f" else if (i <= 4)
    paste0(" f", paste0(rep("'", i), collapse = "")) else paste0(" f^(", i, ")"))
  frac <- paste0(sprintf("%.4e", resulting.terms), flabs)
  frac <- paste0(frac[which(resulting.terms != 0)[1:2]], collapse = " + ")
  frac <- gsub("*1\\.0000e\\+00", "", frac)
  frac <- paste0(gsub("\\+ -", "- ", gsub("^ +", "", frac)), " + ...")
  # Updating the effective accuracy order in case of zero stencil elements
  ea <- diff(which(abs(resulting.terms) > 1024*zero.tol)[1:2])

  ret <- list(stencil = stencil, weights = weights)
  attr(ret, "remainder.coef") <- resulting.terms[which(resulting.terms != 0)[2]]
  attr(ret, "accuracy.order") <- c(requested = acc.requested, effective = ea)
  attr(ret, "expansion") <- frac

  return(ret)
}
