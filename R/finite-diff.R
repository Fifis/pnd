#' Finite-difference coefficients for arbitrary grids
#'
#' This function computes the coefficients that yield an numerical approximation
#' or arbitrary order to the true \eqn{m^{\textrm{th}}}{m-th} derivative of the
#' given function. Given the derivative order and desired accuracy order,
#' returns a grid (stencil) of equally spaced points and the weights for
#' combining the function values obtained on that stencil. Alternatively, for a
#' given arbitrary stencil \eqn{\{b_i\}_{i=1}^n}{{s[i]}, i = 1, ..., n}, it computes
#' the optimal weights \eqn{\{w_i\}}{{w[i]}} yielding
#' the numerical approximation of the derivative:
#' \deqn{\frac{d^m f}{dx^m} \approx h^{-m} \sum_{i=1}^n w_i f(x + b_i\cdot h)}{d^m/dx^m f(x) ~ sum_i w[i] f(x + b[i]*h)}
#'
#' @param deriv.order Order of the derivative (\eqn{m}{m} in \eqn{\frac{d^m f}{dx^m}}{d^m/dx^m f(x)})
#' @param acc.order Order of accuracy in terms of the step size: for accuracy order \eqn{a}{a},
#'     the approximation error is \eqn{O(f^{(a+1)}(c))}{O(d^o/dx^o f(c))}, where \eqn{x -h \le c \le x+h}{x-h <= c <= x+h} (depends on the higher-order derivatives).
#' @param side Character: `"central"` (symmetrical two-sided, default), `"forward"`, or `"backward"` differences.
#'     Unless the function is computationally prohibitively expensive, two-sided differences are strongly recommended.
#' @param stencil In case the user desires a non-uniform grid or a custom grid,
#'     a vector of points at which the function is to be evaluated. For derivative
#'     order `m`, must contain at least `m+1` points.
#' @param zero.action Character: if `"drop"`, stencil points corresponding to
#'     weights less in absolute value than `zero.tol` will be omitted. If `"round"`,
#'     all stencil points are preserved, and small weights (less in absolute
#'     value than `zero.tol`) are replaced with zeros. Otherwise, do nothing.
#'     E.g. the stencil for d/dx f(x) is (-1, 0, 1) with weights (-0.5, 0, 0.5).
#'     Re-computing f(x + 0*h) = f(x) can be wasteful.
#' @param zero.tol Non-negative positive scalar that determines near-zero
#'     weights.
#'
#' @details
#' The finite-difference coefficients for any given stencil are given as a solution of a linear equation
#' system. There derivation of the system is due to \insertCite{taylor2016finite}{pnd}, although a similar
#' approach is described in \insertCite{fornberg1988generation}{pnd}. This function reproduces the tables
#' from the latter paper exactly.
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
#'
#' @examples
#' fdCoef()  # Simple two-sided derivative
#' fdCoef(2) # Simple two-sided second derivative
#' fdCoef(acc.order = 4)$weights * 12  # Should be (1, -8, 8, -1)
#' # Replicating Table 1 from Fornberg (1988) (cited above)
#' pad9 <- function(x) {l <- length(x); c(a <- rep(0, (9-l)/2), x, a)}
#' f <- function(d, a) pad9(fdCoef(deriv.order = d, acc.order = a,
#'                                 zero.action = "round")$weights)
#' t1.1 <- t(sapply((1:4)*2, \(a) f(d = 1, a)))
#' t1.2 <- t(sapply((1:4)*2, \(a) f(d = 2, a)))
#' t1.3 <- t(sapply((1:3)*2, \(a) f(d = 3, a)))
#' t1.4 <- t(sapply((1:3)*2, \(a) f(d = 4, a)))
#' t1 <- data.frame(OrdDer = rep(1:4, times = c(4, 4, 3, 3)),
#'                  OrdAcc = c((1:4)*2, (1:4)*2, (1:3)*2, (1:3)*2),
#'                  rbind(t1.1, t1.2, t1.3, t1.4))
#' colnames(t1)[3:11] <- as.character(-4:4)
#' print(t1, digits = 4)
#'
#' # Replicating Table 2 (ibid) for a halfway stencil (without 0th derivatives)
#' s <- c(-7, -5, -3, -1, 1, 3, 5, 7)/2
#' pad8 <- function(x) {l <- length(x); c(a <- rep(0, (8-l)/2), x, a)}
#' g <- function(d, a) {
#'   k <- a/2 + max(floor(d/2), 0)
#'   pad8(fdCoef(deriv.order = d, stencil = s[(5-k):(4+k)])$weights)
#' }
#' t2.1 <- t(sapply((1:4)*2, \(a) g(d = 1, a)))
#' t2.2 <- t(sapply((1:3)*2, \(a) g(d = 2, a)))
#' t2.3 <- t(sapply((1:3)*2, \(a) g(d = 3, a)))
#' t2.4 <- t(sapply((1:2)*2, \(a) g(d = 4, a)))
#' t2 <- data.frame(OrdDer = rep(1:4, times = c(4, 3, 3, 2)),
#'                  OrdAcc = c((1:4)*2, (1:3)*2, (1:3)*2, (1:2)*2),
#'                  rbind(t2.1, t2.2, t2.3, t2.4))
#' colnames(t2)[3:10] <- as.character(s)
#' print(t2, digits = 4)
#'
#' # Replicating Table 3 (ibid)
#' pad9r <- function(x) c(x, rep(0, 9-length(x)))
#' h <- function(d, a) pad9r(fdCoef(deriv.order = d, stencil = 0:(d+a-1))$weights)
#' t3.1 <- t(sapply(1:8, \(a) h(d = 1, a)))
#' t3.2 <- t(sapply(1:7, \(a) h(d = 2, a)))
#' t3.3 <- t(sapply(1:6, \(a) h(d = 3, a)))
#' t3.4 <- t(sapply(1:5, \(a) h(d = 4, a)))
#' t3 <- data.frame(OrdDer = rep(1:4, times = 8:5),
#'                  OrdAcc = c(1:8, 1:7, 1:6, 1:5),
#'                  rbind(t3.1, t3.2, t3.3, t3.4))
#' colnames(t3)[3:11] <- as.character(0:8)
#' print(t3, digits = 4)
#'
#' # Using an custom stencil for the first derivative: x-2h and x+h
#' fdCoef(stencil = c(-2, 1))
#' @export
fdCoef <- function(deriv.order = 1, side = c("central", "forward", "backward"),
                   acc.order = 2, stencil = NULL,
                   zero.action = c("drop", "round", "none"), zero.tol = (10^deriv.order)*.Machine$double.eps) {
  if (deriv.order > 4) {
    warning(paste0("You are trying to compute the ", deriv.order, "-th derivative -- it can be numerically unstable: the order is too high. Proceed with caution."))
    zero.tol <- 1e-10
  }
  side <- side[1]
  zero.action <- zero.action[1]
  if (length(acc.order) != 1) stop("The 'acc.order' argument must have length 1.")
  if (acc.order < 1 | (side == "central" & acc.order/2 != round(acc.order/2)) | (acc.order != round(acc.order)))
    stop("The order of accuracy must be a positive integer. For 2-sided derivatives, it must be even.")
  l.end <- switch(side,
                  central = acc.order/2 + floor((deriv.order-1)/2),
                  backward = acc.order + deriv.order,
                  forward = acc.order + deriv.order)
  if (is.null(stencil)) stencil <- switch(side,
                                          central = (-l.end):l.end,
                                          backward = (-l.end):0,
                                          forward = 0:l.end)
  l <- length(stencil)
  if (l < deriv.order+1) stop("To compute the d-th derivative, at least d+1 stencil points are required.")
  A <- t(sapply(1:l, function(i) stencil^(i-1)))
  b <- numeric(l); b[1+deriv.order] <- factorial(deriv.order)
  weights <- solve(A, b)
  if (zero.action != "none") {
    zw <- abs(weights) < zero.tol
    if (zero.action == "drop") {
      stencil <- stencil[!zw]
      weights <- weights[!zw]
    } else if (zero.action == "round") {
      weights[zw] <- 0
    } else stop("'zero.action' must be 'drop', 'round', or 'none'")
  }
  rs <- round(stencil, 2) # For names
  names(weights) <- paste0("x", ifelse(rs < 0, "-", "+"), abs(rs), "h")
  names(weights)[rs == 0] <- "x"
  return(list(stencil = stencil, weights = weights))
}


#' Parallelised gradient computation
#'
#' Computes a two- or one-sided numerical derivative that approximates the gradient | Jacobian using the indicated number of cores for maximum efficiency.
#'
#' @param func A function that returns a numeric scalar or a vector. If the function is vector-valued, the, the result is the Jacobian.
#' @param x A point at which the gradient or Jacobian needs to be estimated.
#' @param h The numerical difference step size. Too large = the slope of the secant is a bad estimator of the gradient, too small = ill conditioning (0/0).
#' @param acc.order Desired order of accuracy. The error is usually O(h^acc.order). To achieve this order of accuracy, the function needs to be evaluated \code{acc.order*length(x)} times.
#' @param side Passed to \code{fdCoef()}. Centred or one-sided differences. Unless the function is computationally prohibitively expensive, two-sided differences are strongly recommended.
#' @param parallel If TRUE, estimates the gradient via finite differences where the function is evaluated in parallel.
#' @param cores Number of forked processes.
#' @param cluster A cluster on which the computations are done.
#' @param load.balance If TRUE, disables pre-scheduling for \code{mclapply} or enables load balancing via \code{parLapplyLB}.
#' @param ... Passed to `func`.
#'
#' Note that for one-sided problems, the step size that make the formula error
#' equal to the truncation error is of the order Mach.eps^(1/2) and for two-sided, Mach.eps^(1/3).
#' However, the optimal step size depends on the value of the higher-order derivatives
#' that is not available in general (or required extra computation that is, in turn, prone to numerical error).
#'
#' @return If \code{func} returns a scalar, a vector of the same length as \code{x}.
#' If \code{func} returns a vector, then, a matrix of dimensions \code{length(f(x)) length(x)}
#'
#' @examples
#' \dontrun{
#' slowFunScalar <- function(x) {Sys.sleep(0.04); print(x, digits = 12); sum(sin(x))}
#' slowFunVector <- function(x) {Sys.sleep(0.04); print(x, digits = 12); c(sum(sin(x)), sum(exp(x)))}
#' true.g <- cos(1:4) # Analytical gradient
#' true.j <- rbind(cos(1:4), exp(1:4)) # Analytical Jacobian
#' system.time(g.slow <- numDeriv::grad(slowFunScalar, x = 1:4) - true.g)
#' system.time(j.slow <- numDeriv::jacobian(slowFunVector, x = 1:4) - true.j)
#' system.time(g.fast <- Grad(slowFunScalar, x = 1:4,
#'                                    parallel = TRUE, cores = 4) - true.g)
#' system.time(j.fast <- Grad(slowFunVector, x = 1:4,
#'                                    parallel = TRUE, cores = 4) - true.j)
#' system.time(j.fast4 <- Grad(slowFunVector, x = 1:4, acc.order = 4,
#'                                     parallel = TRUE, cores = 4) - true.j)
#' rownames(j.slow) <- c("numDeriv.jacobian", "")
#' rownames(j.fast) <- c("fast.jacobian.order2", "")
#' rownames(j.fast4) <- c("fast.jacobian.order4", "")
#' # Discrepancy
#' rbind(numDeriv.grad = g.slow, fast.grad = g.fast, j.slow, j.fast, j.fast4)
#' # The order-4 derivative is more accurate for functions with large high-order derivatives
#'}
#'
#' @export
Grad <- function(func, x, side = c("central", "forward", "backward"),
                 acc.order = if (side[1] == "central") 2 else 1,
                 h = .Machine$double.eps^(1/(1+acc.order)),
                 parallel = c("auto", "fork", "PSOCK", "none"), cores = 2, cluster = NULL, load.balance = TRUE,
                 ...) {
  parallel <- parallel[1]
  parallel <- parallel != "none"
  side <- side[1]
  n <- length(x)
  s <- fdCoef(acc.order = acc.order, side = side)

  # Parallelising the task in the most efficient as possible: we need as many evaluations
  # as length(x) * acc.order
  dx <- s$stencil*h
  par.grid <- expand.grid(arg = 1:n, grid = 1:length(dx))
  par.grid$step <- dx[par.grid$grid]
  par.grid$weights <- s$weights[par.grid$grid]
  xx <- lapply(1:nrow(par.grid), function(i) {
    newx <- x
    newx[par.grid$arg[i]] <- newx[par.grid$arg[i]] + par.grid$step[i]
    newx
  })
  ff <- if (parallel & cores > 1) {
    parallel::mclapply(X = xx, FUN = \(x) func(x, ...), mc.cores = cores, mc.preschedule = !load.balance)
  } else {
    lapply(xx, \(x) func(x, ...))
  }
  # else if (parallel & !use.mclapply & !is.null(cluster)) {
  #if (load.balance) parallel::parLapplyLB(cl = cluster, X = xx, fun = \(x) func(x, ...)) else
  #    parallel::parLapply(cl = cluster, X = xx, fun = \(x) func(x, ...))
  # The output can be vector-valued, which is why we work carefully with rows
  ff <- do.call(rbind, ff)
  ffw <- ff * par.grid$weights
  ffl <- split(as.data.frame(ffw), f = par.grid$arg)
  if (ncol(ff) == 1) ffl <- lapply(ffl, as.matrix) # The dimensions must be preserved
  ffs <- lapply(ffl, colSums) # A list of vectors or scalars
  jac <- unname(as.matrix(do.call(cbind, ffs)))
  if (nrow(jac) == 1) {
    jac <- as.numeric(jac)
    if (!is.null(names(x))) names(jac) <- names(x)
  } else {
    if (!is.null(names(x))) colnames(jac) <- names(x)
  }
  return(jac / h)
}

# rm(list = ls())
# fScal <- \(x) {Sys.sleep(0.04); print(x, 12); sum(sin(x))}
# fVect <- \(x) {Sys.sleep(0.04); print(x, 12); c(sum(sin(x)), sum(exp(x)))}
# true.g <- cos(1:4) # Analytical gradient
# true.j <- rbind(cos(1:4), exp(1:4)) # Analytical Jacobian
# system.time(discr.g.slow <- numDeriv::grad(fScal, x = 1:4) - true.g)
# system.time(discr.j.slow <- numDeriv::jacobian(fVect, x = 1:4) - true.j)
# system.time(discr.g.fast <- Grad(fScal, x = 1:4, cores = 4) - true.g)
# system.time(discr.j.fast <- Jacobian(fVect, x = 1:4, cores = 4) - true.g)
