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
#' A valid truncation branch requires 3 consecutive centred slopes within 10%
#' of \code{acc.order}, without bridging missing estimates. The fit uses all finite
#' estimates up to the end of the longest such run.
#' A failed fit or a kink with observations on only one side invokes the fallback.
#'
#' The accuracy of the estimated higher-order derivatives starts at 2;
#' only zero or non-finite estimates are retried at accuracy 4.
#' Auxiliary weighted sums within their arithmetic rounding bound are treated as
#' unresolved rather than positive observations in the fit.
#' All stencils reuse the evaluated grid; no extra steps are evaluated for retries.
#'
#' Interpolation uses the existing derivative estimates and their leading dependence
#' on \code{h^acc.order}, avoiding further function evaluations. For continuous
#' functions, \eqn{f(x+\delta)\approx f(x)} for small \eqn{\delta}.
#'
#' The estimated step must lie in the search interval. Interpolation requires two
#' finite neighbouring derivative estimates; otherwise the nearest usable grid step
#' is returned with its derivative estimate. Steps with infinite or coincident
#' stencil arguments are excluded. The function stops if no usable grid step remains.
#' These checks reuse the evaluated grid and make no further function calls.
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
#'         \item \code{3} – The fitted kink has observations on only one side; a fallback is returned.
#'         \item \code{4} – Step trimmed to 0.1|x| when |x| is not tiny and within range.
#'         \item \code{5} – The fit or selected step failed a safeguard; a fallback is returned.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the step and argument grids,
#'       function values on those grids, estimated derivative values, estimated error values,
#'       and predicted model-based errors. The \code{fhigher} column contains the higher-order
#'       auxiliary derivative; \code{aprime} records its selected accuracy order.
#'       \code{fit} contains the initial and fitted log-kink coordinates (\code{gamma}, \code{beta}),
#'       unscaled initial residual \code{mad}, fixed knee \code{rho} (at least one binary-logarithm unit), fitting \code{subset},
#'       mean pseudo-Huber \code{loss}, optimiser \code{convergence} and whether the fit was \code{accepted}; a \code{NULL}
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
    range <- h0 * c(2^-36, 2^24)
  }
  if (length(range) != 2 || any(!is.finite(range)) || any(range <= 0))
    stop("The range must be a positive vector of length 2 with finite endpoints.")
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
                            preschedule = preschedule, check.zero = TRUE, ...)
  fgrid <- estimates$f
  cds <- estimates$deriv
  fhigher <- estimates$higher
  K <- estimates$K
  etrunc <- estimates$trunc
  eround <- estimates$round
  fallback <- NULL

  log2h <- log2(hgrid)
  log2e <- suppressWarnings(log2(etrunc))
  log2e[!is.finite(log2e)] <- NA_real_

  # The two straight branches meet at theta[1], with height theta[2]
  checkFit <- function(theta, x, a, d)
    theta[2] + ifelse(x < theta[1], -d * (x-theta[1]), a * (x-theta[1]))

  i.finite <- which(is.finite(log2e))
  i.trunc.valid <- i.round <- integer(0)
  if (length(i.finite) < 3L) {
    exitcode <- 2L
  } else {
    # A centred slope needs three consecutive observations; do not bridge missing estimates
    middle <- seq.int(2L, n-1L)
    complete <- is.finite(log2e[middle-1L]) & is.finite(log2e[middle]) & is.finite(log2e[middle+1L])
    local.slopes <- rep(NA_real_, n)
    local.slopes[middle[complete]] <-
      (log2e[middle[complete]+1L] - log2e[middle[complete]-1L]) /
      (log2h[middle[complete]+1L] - log2h[middle[complete]-1L])
    good <- which(abs(local.slopes/acc.order - 1) < 0.1)
    if (length(good)) {
      runs <- splitRuns(good)
      i.trunc.valid <- runs[[which.max(lengths(runs))]]
    }

    if (length(i.trunc.valid) < 3L) {
      warning("Fewer than 3 consecutive observations in the truncation branch! Resorting to a safe option.")
      exitcode <- 1L
    } else {
      # Include the smaller steps to identify the rounding branch as well as the truncation branch
      sbst <- i.finite[i.finite <= max(i.trunc.valid)]
      i.round <- i.finite[i.finite < min(i.trunc.valid)]

      # The median log-step among the 5 lowest errors prevents fragmented local optima
      # as initial values; this is a purely empirical improvement
      i.low <- utils::head(sbst[order(log2e[sbst])], 5)
      theta0 <- c(stats::median(log2h[i.low]), min(log2e[sbst]))
      if (theta0[1] > log2h[max(i.trunc.valid)])
        theta0 <- c(log2h[min(i.trunc.valid)], log2e[min(i.trunc.valid)])

      resid0 <- log2e[sbst] - checkFit(theta0, log2h[sbst], acc.order, deriv.order)
      rmad <- stats::mad(resid0, constant = 1)
      # Keep the knee fixed; use unity (one power of two) instead for near-zero MADs
      rho <- if (is.finite(rmad)) max(1, rmad) else 1
      fobj <- function(theta) {
        resid <- log2e[sbst] - checkFit(theta, log2h[sbst], acc.order, deriv.order)
        # Equivalent pseudo-Huber loss, avoiding cancellation and division by a tiny rho
        mean(rho * resid^2 / (sqrt(rho^2 + resid^2) + rho))
      }
      opt <- tryCatch(stats::optim(par = theta0, fn = fobj, method = "BFGS",
                                  control = list(reltol = .Machine$double.eps^(1/3))),
                      error = function(e) NULL)
      fit <- list(initial = stats::setNames(theta0, c("gamma", "beta")),
                  par = if (is.null(opt)) c(gamma = NA_real_, beta = NA_real_)
                        else stats::setNames(opt$par, c("gamma", "beta")),
                  mad = rmad, rho = rho, subset = sbst,
                  loss = if (is.null(opt)) NA_real_ else opt$value,
                  convergence = if (is.null(opt)) NA_integer_ else opt$convergence,
                  accepted = FALSE)

      if (is.null(opt) || opt$convergence != 0L ||
          any(!is.finite(opt$par)) || !is.finite(opt$value)) {
        exitcode <- 5L
      } else if (!any(log2h[sbst] < opt$par[1]) || !any(log2h[sbst] > opt$par[1])) {
        # A one-sided fit cannot identify the balance between the two error branches
        exitcode <- 3L
      } else {
        # The two-branch model has T/R = d/a at its minimum, giving h* = hcross*(d/a)^(1/(d+a))
        # TODO: calibrate this heuristic for the auxiliary/target rounding-bound ratio lambda
        # (11/16 for the 6-point stencil, 1/2 for 4 points at d=1, a=2)
        hopt <- 2^opt$par[1] * (deriv.order/acc.order)^(1/(deriv.order + acc.order))
        if (!is.finite(hopt) || hopt < range[1] || hopt > range[2]) {
          exitcode <- 5L
        } else {
          fit$accepted <- TRUE
          hat.check <- 2^checkFit(opt$par, x = log2h, a = acc.order, d = deriv.order)
        }
      }
    }
  }

  if (exitcode %in% c(1L, 2L, 3L, 5L)) {
    fallback <- fallbackM(hgrid, fhigher, K, estimates$round.coef, eround,
                          x, deriv.order, acc.order)
    hopt <- fallback$par
  }

  # Keep the existing local-step cap only when it belongs to the permitted interval
  upper <- range[2]
  cap <- 0.1 * abs(x)
  if (abs(x) > 4.71216091538e-7 && cap >= range[1] && cap <= upper) {
    upper <- cap
    if (is.finite(hopt) && hopt > cap) {
      exitcode <- 4L
      hopt <- cap
    }
  }

  # A finite cached stencil above hopt encloses its target nodes on an interval domain
  # Reject coincident or infinite arguments even if FUN happens to return a finite number there
  target.x <- x + outer(hgrid, estimates$stencil)
  finite.args <- rowSums(!is.finite(target.x)) == 0L
  # Nodes are ordered, so any arguments that coincide must be adjacent
  distinct.args <- rowSums(target.x[, -1L, drop = FALSE] ==
                             target.x[, -ncol(target.x), drop = FALSE]) == 0L
  valid.target <- which(is.finite(cds) & finite.args & distinct.args)
  available <- valid.target[hgrid[valid.target] >= range[1] & hgrid[valid.target] <= upper]
  if (!length(available)) stop("Could not compute a finite derivative on the permitted grid.")
  neighbours <- if (is.finite(hopt) && hopt > 0)
    utils::head(order(abs(log2h - log2(hopt))), 2L) else integer(0)
  exact <- which(hgrid == hopt)
  supported <- is.finite(hopt) && hopt >= range[1] && hopt <= upper &&
    if (length(exact)) exact[1L] %in% available else
      all(neighbours %in% valid.target) && any(hgrid[valid.target] >= hopt) &&
      all(is.finite(x + estimates$stencil*hopt)) &&
      !anyDuplicated(x + estimates$stencil*hopt)
  if (!supported) {
    # Return an evaluated step when interpolation cannot be supported, keeping par and value consistent
    distance <- if (is.finite(hopt) && hopt > 0) abs(log2h[available] - log2(hopt)) else rep(0, length(available))
    hopt <- hgrid[available[which.min(distance)]]
    exitcode <- 5L
  }
  log2hopt <- log2(hopt)

  cd <- interpolateM(hgrid, cds, hopt, acc.order)
  if (!is.finite(cd)) stop("Could not compute a finite derivative on the grid.")
  if (!is.null(fit) && isTRUE(fit$accepted)) {
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
                "the fitted kink has observations on only one side, using a rounding-error fallback",  # 3
                "step size too large relative to x, using |x|/10 instead",  # 4
                "fit or selected step failed a safeguard, using a fallback"  # 5

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
