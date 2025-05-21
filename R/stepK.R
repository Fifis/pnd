#' Kink-based step selection
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
#' plot(step.K(sin, 1))
#' step.K(exp, 1, range = c(1e-12, 1e-0))
#' step.K(atan, 1)
#'
#' # Edge case 1: function symmetric around x0, zero truncation error
#' step.K(sin, pi/2)
#' step.K(sin, pi/2, shrink.factor = 0.8)
#'
#' # Edge case 1: the truncation error is always zero and f(x0) = 0
#' suppressWarnings(step.K(function(x) x^2, 0))
#' # Edge case 2: the truncation error is always zero
#' step.K(function(x) x^2, 1)
#' step.K(function(x) x^4, 0)
#' step.K(function(x) x^4, 0.1)
#' step.K(function(x) x^6 - x^4, 0.1)
#' step.K(atan, 3/4)
#' step.K(exp, 2, range = c(1e-16, 1e+1))
step.K <- function(FUN, x, h0 = NULL, deriv.order = 1, acc.order = 2,
                   range = NULL, shrink.factor = 0.5,
                   max.rel.error = .Machine$double.eps^(7/8),
                   cores = 1, preschedule = getOption("pnd.preschedule", TRUE), cl = NULL, ...) {
  if (length(x) != 1) stop(paste0("The step-size selection can handle only univariate inputs. ",
                                  "For 'x' longer than 1, use 'gradstep'."))
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
    warning("The initial range [", toString(printE(range.old)), "] was extended to [",
            toString(printE(range)), "] to ensure a large-enough search space.")
  }

  exitcode <- 0L

  # Creating a sequence of step sizes for evaluation
  g <- gridM(x = x, range = range, shrink.factor = shrink.factor)
  hgrid <- g$h
  xgrid <- g$x
  n <- length(hgrid)
  hat.check <- rep(NA, length(hgrid))

  fgrid <- getValsM(FUN = FUN, x = xgrid, cores = cores, cl = cl, preschedule = preschedule, ...)
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
  cds <- vector("list", 3)
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
  i.trunc.valid <- i.round <- NULL
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

      # Checking monotonicity after this addition to tackle edge cases like
      # slopes = c(-1, 2, -1, -1, 2, -1) due to merging
      # there are at least 3 points for sure, so the middle ones need to be taken out
      nonmonot.i <- which(local.slopes[i.trunc.valid[2:(length(i.trunc.valid)-1)]] < 0) + 1
      if (length(nonmonot.i) > 0) i.trunc.valid <- i.trunc.valid[-nonmonot.i]

      # Another edge case: if the slope of etrunc is always acc.order, then, there is no check
      if (0 %in% i.trunc.valid) {
        only.right.branch <- FALSE
        i.trunc.valid <- i.trunc.valid[i.trunc.valid > 0]
        exitcode <- 3L
      }

      # If there are more than 1 run, pick the longest one
      # However, if it is shorter than 5 (2 points were added artificially; we want a length-3 run), discard it
      itv.runs <- splitRuns(i.trunc.valid)
      itv.len <- sapply(itv.runs, length)
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
      warning("Fewer than 5 observations in the truncation branch! Resorting to a safe option.")
      # At the optimum, the rounding error should be approx. |(e^(2/3) f^(2/3) f'''^(1/3))/(2^(2/3) 3^(1/3))|
      # However, if f''' ~ 0, we use only the rounding branch and compare it to the optimal value
      # if f''' ~ 1
      exitcode <- 1L
    }
  } else {
    zero.trunc <- TRUE
    exitcode <- 2L
  }

  if (exitcode %in% c(1L, 2L)) {
    complete.rows <- apply(fgrid, 1, function(x) all(is.finite(x)))
    f0 <- abs(max(abs(fgrid[which(complete.rows)[1], ])))
    if (!is.finite(f0)) stop(paste0("Could not compute the function value at ", x, "."))
    if (f0 < .Machine$double.eps) f0 <- .Machine$double.eps
    expected.eround <- (.Machine$double.eps^2 * f0^2 / 12)^(1/3)
    # TODO: generalise
    # Two cases possible: f(x) = x^2 at x = 0 has eround increasing,
    # which implies that a large fail-safe step can be chosen
    f.er <- is.finite(eround)
    mean.sign <- if (sum(f.er) > 1) mean(sign(diff(eround[f.er]))) else 1
    # If the rounding error is growing, prefer a slightly larger step size
    if (mean.sign > 0.5) {
      hopt <- 128*stepx(x = x, deriv.order = deriv.order, acc.order = acc.order,
                        zero.tol = .Machine$double.eps^(1/3))
    } else {
      # Otherwise, for an eround that does not grow convincingly,
      # find the step that approximates the theoretical expected value
      hopt  <- hgrid[which.min(abs(eround - expected.eround))]
    }
    log2hopt <- log2(hopt)
    # TODO: generalise; do the theory!
  } else if (exitcode == 3L) {  # If the right branch is growing, choose a smaller step size
    hopt <- hgrid[round(stats::quantile(i.trunc.valid, 0.25))]
    log2hopt <- log2(hopt)
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

  msg <- switch(exitcode + 1L,
                "target error ratio reached within tolerance",
                "truncation branch slope is close to zero, relying on the expected rounding error",
                "truncation error is near-zero, relying on the expected rounding error",
                "there is no left branch of the combined error, fitting impossible"

  )

  diag.list <- list(h = hgrid, x = xgrid, f = fgrid, deriv = cbind(f = cds[[1]], fhigher = cds[[2]]),
                    est.error = cbind(trunc = etrunc, round = eround, check = hat.check))

  ret <- list(par = hopt, value = cd, counts = n, exitcode = exitcode, message = msg,
              abs.error = c(trunc = et, round = er), method = "Kostyrka", iterations = diag.list,
              args = list(deriv.order = deriv.order, acc.order = acc.order,
                          range = range, shrink.factor = shrink.factor, max.rel.error = max.rel.error,
                          i.good = i.trunc.valid, i.okay = NULL, i.increasing = NULL, i.round = i.round))
  class(ret) <- "stepsize"
  return(ret)
}


plotK <- function(x, ...) {
  cols <- c("#7e1fde", "#328d2d") # , "#d58726", "#ca203a")
  it <- x$iterations
  et <- it$est.error[, "trunc"]
  er <- it$est.error[, "round"]
  et[et == 0] <- .Machine$double.eps
  er[er == 0] <- .Machine$double.eps
  if (!any(is.finite(et) & (et != 0))) {
    warning("The truncation error is either exactly 0 or NA. Nothing to plot.")
    return(invisible(x))
  }
  h <- it$h
  xl <- range(h)
  xl <- sqrt(range(h) * xl)
  yl <- range(et[is.finite(et)]) * c(0.5, 2)
  plot(h[et > 0], et[et > 0],
       log = "xy", bty = "n", ylim = yl,
       ylab = "Estimated abs. error in df/dx", xlab = "Step size",
       main = "Kostyrka step-size selection", ...)
  # graphics::mtext(paste0("dotted line: rounding error assuming max. rel. err. < ", printE(epsilon, 1)), cex = 0.8, line = 0.5)
  if (length(x$args$i.round) > 0)
    graphics::points(h[x$args$i.round], et[x$args$i.round], pch = 16, col = cols[1], cex = 0.9)
  # if (length(x$args$i.increasing) > 0)  # This is always NULL -- for compatibility with plotM
  #   graphics::points(h[x$args$i.increasing], et[x$args$i.increasing], pch = 16, col = cols[4], cex = 0.9)
  # if (length(x$args$i.ok) > 0)  # This is always NULL -- for compatibility with plotM
  #   graphics::points(h[x$args$i.ok], et[x$args$i.ok], pch = 16, col = cols[3], cex = 0.9)
  if (length(x$args$i.good) > 0)
    graphics::points(h[x$args$i.good], et[x$args$i.good], pch = 16, col = cols[2], cex = 0.9)
  i.bad <- setdiff(seq_along(h), c(x$args$i.good, x$args$i.ok, x$args$i.increasing, x$args$i.round))
  if (length(i.bad) > 0)
    graphics::points(h[i.bad], et[i.bad], pch = 4, col = "#00000088", cex = 0.9)

  # graphics::lines(h, er, col = "#000000", lty = 3)

  graphics::lines(h, it$est.error[, "check"], col = "#FFFFFFBB", lwd = 3)
  graphics::lines(h, it$est.error[, "check"])

  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Rounding", "Trunc. good", "Invalid"),
                   pch = c(16, 16, 4), col = c(cols, "#000000"),
                   box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(x)
}
