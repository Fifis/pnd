#' Kink-based step selection
#'
#' Optimal step-size search by fitting the rounding and truncation branches of
#' observed error estimates and locating the transition between them.
#'
#' @inheritParams step.M
#' @inheritParams GenD
#'
#' @details
#' This function estimates the optimal step size by fitting a V-shaped curve to error estimates on log--log axes,
#' i.e. finding a regression kink, using
#' unconstrained BFGS optimisation and a mean pseudo-Huber penalty, which limits
#' the influence of irregular estimates. The rounding and truncation branches have
#' slopes \code{-deriv.order} and \code{acc.order}; for two-point central differences,
#' these are -1 and 2, respectively. The step at the fitted kink is multiplied by
#' \eqn{(d/a)^{1/(d+a)}}, where \eqn{d} is \code{deriv.order} and \eqn{a} is
#' \code{acc.order}, to estimate the minimum of the combined error.
#' If too few estimates follow the expected truncation slope, a rounding-error
#' fallback is used instead.
#'
#' The accuracy of the estimated higher-order derivatives starts at 2;
#' only zero or non-finite estimates are retried at accuracy 4.
#' All stencils reuse the evaluated grid; no extra steps are evaluated for retries.
#'
#' Interpolation uses the existing derivative estimates and their leading dependence
#' on \code{h^acc.order}, avoiding further function evaluations. For continuous
#' functions, \eqn{f(x+\delta)\approx f(x)} for small \eqn{\delta}.
#'
#' If only one neighbouring grid point has a finite derivative estimate, it is reused.
#' If neither neighbour has a finite estimate, the finite derivative estimate closest
#' to the middle of the step-size grid is used. The function stops with an error if
#' no finite derivative estimate is available.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical derivative (using central differences).
#'     \item \code{counts} – the number of paired grid points, with one additional central evaluation for even derivative orders.
#'     \item \code{abs.error} – separate fitted error components, or interpolated grid estimates on fallback;
#'       the estimated rounding component describes noise in the truncation-error estimates.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – A V-shaped curve was fitted without errors or sanity-check failures.
#'         \item \code{1} – The truncation branch is ill-behaved; a fail-safe value is returned.
#'         \item \code{2} – Too few positive finite estimates; a fail-safe value is returned.
#'         \item \code{3} – There is no left branch of the V shape; a fail-safe value is returned.
#'         \item \code{4} – Step trimmed to 0.1|x| when |x| is not tiny and within range.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the step and argument grids,
#'       function values on those grids, estimated derivative values, estimated error values,
#'       and predicted model-based errors. The \code{fhigher} column contains the higher-order
#'       auxiliary derivative; \code{aprime} records its selected accuracy order.
#'       \code{fit} contains the initial and fitted log-kink coordinates (\code{gamma}, \code{beta}),
#'       unscaled initial residual \code{mad}, fixed knee \code{rho} (at least one binary-logarithm unit), fitting \code{subset},
#'       mean pseudo-Huber \code{loss} and optimiser \code{convergence}; a \code{NULL}
#'       is returned if no fit was attempted.
#'       \code{fallback} records the quantities used when the truncation branch cannot be found.
#'       Its \code{target} is the rounding error predicted at the minimum of the combined-error model.
#'       It is \code{NA} when the estimates do not identify a finite minimum.
#'       The calculation uses the median absolute higher-derivative estimate and the function values
#'       at the smallest step for which all points in the difference formula are available.
#'       \code{mean.sign} indicates whether rounding estimates tend to rise or fall as steps increase.
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
#' # Edge case 2: the truncation error is always zero and f(x0) = 0
#' suppressWarnings(step.K(function(x) x^2, 0))
#' # Other edge cases
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
  if (length(x) != 1) stop(paste0("Direct step-size selection can handle only univariate inputs. ",
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

  cl <- getParallelCluster(cores, cl)
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
  g <- gridM(x = x, range = range, ratio = inv.sf)
  hgrid <- g$h
  xgrid <- g$x
  n <- length(hgrid)
  hat.check <- rep(NA, length(hgrid))
  fit <- NULL

  estimates <- differencesM(FUN, x, g, deriv.order, acc.order, inv.sf, retry = TRUE,
                            max.rel.error = max.rel.error, cores = cores, cl = cl,
                            preschedule = preschedule, ...)
  fgrid <- estimates$f
  cds <- estimates$deriv
  fhigher <- estimates$higher
  K <- estimates$K
  etrunc <- estimates$trunc
  eround <- estimates$round
  fallback <- NULL

  # A value to compare to the threshold for Zero Truncation Error
  good.te <- etrunc != 0 & is.finite(etrunc)
  med.te <- if (any(good.te)) stats::median(etrunc[good.te]) else .Machine$double.eps

  # Initial value:
  log2h <- log2(hgrid)
  log2e <- suppressWarnings(log2(etrunc))
  log2e[is.infinite(log2e)] <- NA  # Preserving NaNs just in case
  # plot(log2h, log2e)

  # Check-like function of two parameters: horizontal and vertical location
  checkFit <- function (theta, x, a, d) (theta[2]-d*(x-theta[1]))*(x<theta[1]) + (theta[2]+a*(x-theta[1]))*(x>=theta[1])

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
        i.round <- seq_len(min(i.trunc.valid) - 1L)
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
      sbst <- intersect(sort(unique(c(i.trunc.valid, i.round))), which(is.finite(log2e)))

      # The median log-step among the 5 lowest errors prevents fragmented local optima
      # as initial values; this is a purely empirical improvement
      i.low <- utils::head(sbst[order(log2e[sbst])], 5)
      theta0 <- c(stats::median(log2h[i.low]), min(log2e[sbst]))
      if (theta0[1] > log2h[max(i.trunc.valid)]) theta0 <- c(log2h[min(i.trunc.valid)], log2e[min(i.trunc.valid)])
      # If the rounding error was not available, use a rough fall-back
      if (is.na(theta0[2])) theta0[2] <- stats::quantile(log2e, 0.25, na.rm = TRUE)

      resid0 <- log2e[sbst] - checkFit(theta0, log2h[sbst], acc.order, deriv.order)
      rmad <- stats::mad(resid0, constant = 1)
      # Keep the knee fixed; use unity (one power of two) instead for near-zero MADs
      rho <- if (is.finite(rmad)) max(1, rmad) else 1
      fobj <- function(theta) {
        resid <- log2e[sbst] - checkFit(theta, log2h[sbst], acc.order, deriv.order)
        # Equivalent pseudo-Huber loss, avoiding cancellation and division by a tiny rho
        mean(rho * resid^2 / (sqrt(rho^2 + resid^2) + rho))
      }
      opt <- stats::optim(par = theta0, fn = fobj, method = "BFGS", control = list(reltol = .Machine$double.eps^(1/3)))
      fit <- list(initial = stats::setNames(theta0, c("gamma", "beta")),
                  par = stats::setNames(opt$par, c("gamma", "beta")),
                  mad = rmad, rho = rho, subset = sbst, loss = opt$value, convergence = opt$convergence)

      # The two-branch model has T/R = d/a at its minimum, giving h* = hcross*(d/a)^(1/(d+a))
      # TODO: calibrate this heuristic for the auxiliary/target rounding-bound ratio lambda
      # (11/16 for the 6-point stencil, 1/2 for 4 points at d=1, a=2)
      hopt0   <- if (!zero.trunc) 2^opt$par[1] else 0.001 * max(1, abs(x))
      hopt <- hopt0 * (deriv.order/acc.order)^(1/(deriv.order + acc.order))
      hat.check <- 2^checkFit(opt$par, x = log2h, a = acc.order, d = deriv.order)
    } else {  # Zero truncation --> choosing the smallest error that yields a non-zero total-error estimate
      hopt <- hgrid[min(i.etrunc.finite)]
    }
    log2hopt <- log2(hopt)
    # points(log2h[sbst], log2e[sbst], pch = 16)
    # lines(log2h, log2(hat.check), lty = 2)

    if (length(i.trunc.valid) < 5) {
      warning("Fewer than 5 observations in the truncation branch! Resorting to a safe option.")
      # An unresolved branch uses the general model target or an explicit rough fallback
      exitcode <- 1L
    }
  } else {
    zero.trunc <- TRUE
    exitcode <- 2L
  }

  if (exitcode %in% c(1L, 2L)) {
    fallback <- fallbackM(hgrid, fhigher, K, estimates$round.coef, eround,
                          x, deriv.order, acc.order)
    hopt <- fallback$par
    log2hopt <- log2(hopt)
  } else if (exitcode == 3L) {  # If the right branch is growing, choose a smaller step size
    hopt <- hgrid[round(stats::quantile(i.trunc.valid, 0.25))]
    log2hopt <- log2(hopt)
  }

  # Locate neighbouring steps before applying the existing size safeguard
  i.closest <- which(rank(abs(log2(hgrid) - log2hopt), ties.method = "first") <= 2)
  h.closest <- hgrid[i.closest]

  # Final step size check: not too large
  if (max(h.closest) > 0.1*abs(x) && abs(x) > 4.71216091538e-7) {
    exitcode <- 4L
    hopt <- 0.1*abs(x)
    log2hopt <- log2(hopt)
    i.closest <- which(rank(abs(log2(hgrid) - log2hopt), ties.method = "first") <= 2)
    h.closest <- hgrid[i.closest]
  }

  cd <- interpolateM(hgrid, cds, hopt, acc.order)
  if (!is.finite(cd)) stop("Could not compute a finite derivative on the grid.")
  if (!is.null(fit) && !(exitcode %in% c(1L, 2L))) {
    # These are the separate empirical branches, without the deferred lambda calibration
    et <- 2^(opt$par[2] + acc.order * (log2hopt - opt$par[1]))
    er <- 2^(opt$par[2] - deriv.order * (log2hopt - opt$par[1]))
  } else {
    et <- interpolateM(hgrid, etrunc, hopt, acc.order)
    er <- interpolateM(hgrid, eround, hopt, -deriv.order)
  }

  msg <- switch(exitcode + 1L,
                "target error ratio reached within tolerance",  # 0
                "truncation branch unresolved, using a rounding-error fallback",  # 1
                "too few positive finite truncation estimates, using a rounding-error fallback",  # 2
                "there is no left branch of the combined error, fitting impossible",  # 3
                "step size too large relative to x, using |x|/10 instead"  # 4

  )

  diag.list <- list(h = hgrid, x = xgrid, f = fgrid,
                   deriv = cbind(f = cds, fhigher = fhigher), aprime = estimates$aprime,
                   est.error = cbind(trunc = etrunc, round = eround, check = hat.check),
                   fit = fit, fallback = fallback, f0 = estimates$f0)

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
  i.missing <- which(!is.finite(et) | et <= 0)
  et[i.missing] <- NA_real_
  if (!any(is.finite(et))) {
    warning("The truncation error is either exactly 0 or NA. Nothing to plot.")
    return(invisible(x))
  }
  h <- it$h
  yl <- range(et[is.finite(et)]) * c(0.5, 2)
  plot(h[is.finite(et)], et[is.finite(et)],
       log = "xy", bty = "n", xlim = range(h), ylim = yl,
       ylab = "Estimated truncation error", xlab = "Step size",
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
  graphics::points(h[i.missing], rep(min(et, na.rm = TRUE), length(i.missing)),
                   pch = 4, col = "#00000088", cex = 0.9)

  # graphics::lines(h, er, col = "#000000", lty = 3)

  graphics::lines(h, it$est.error[, "check"], col = "#FFFFFFBB", lwd = 3)
  graphics::lines(h, it$est.error[, "check"])

  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Rounding", "Trunc. good", "Invalid"),
                   pch = c(16, 16, 4), col = c(cols, "#000000"),
                   box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(x)
}
