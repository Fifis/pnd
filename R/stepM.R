# Generator for Mathur's and Kostyrka's methods that evaluates the shared function grid
# once so all difference formulae reuse the same values
getValsM <- function(FUN, x, cores, cl, preschedule, ...) {
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = x, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)
  return(fgrid)
}


#' Mathur's AutoDX-like automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @inheritParams GenD
#' @param h0 Numeric scalar: initial step size, defaulting to
#'   \code{2^round(log2(0.001 * max(abs(x), 1)))}.
#' @param range Positive finite vector of length 2 defining the initial search interval.
#'   A narrow interval is extended, with a warning, to cover steps on both sides
#'   of the usual machine-precision estimate. The interval used is saved in \code{args$range}.
#' @param max.rel.error Error bound for the relative function-evaluation error
#'   (\eqn{\frac{\hat f(\hat x) - f(x)}{f(x)}}{(^f(^x) - f(x))/f(x)}). Measures how noisy a function is.
#'   If the function is relying on numerical optimisation routines, consider setting to
#'   \code{sqrt(.Machine$double.eps)}.
#'   If the function has full precision to the last bit, set to \code{.Machine$double.eps/2}.
#' @param shrink.factor A positive scalar less than 1 that is used to create a sequence of
#'   step sizes. The recommended value is 0.5. Change to 0.25 for a faster search. This
#'   number should be a negative power of 2 for the most accurate representation.
#' @param seq.tol Numeric scalar: maximum relative deviation of truncation-error
#'   slopes from \code{acc.order}.
#' @param min.valid.slopes Positive integer: how many points must form a sequence
#'   with the correct slope with relative difference from \code{acc.order} less than \code{seq.tol}.
#'   If \code{shrink.factor} is small (< 0.33), consider reducing this to 4.
#' @param correction Logical: if \code{TRUE}, adjust the selected grid step to achieve the
#'   optimal truncation-to-rounding-error ratio (see Details); otherwise return the selected grid step.
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{mathur2012analytical}{pnd} algorithm. It consists of the following steps.
#'
#' 1. Choose a reasonably large (but not too large) initial step size \eqn{h_0}{h0} and a reduction
#' factor (1/2 for fast, 1/4 for slow functions is a reasonable choice).
#' 2. Estimate truncation errors using a higher-order derivative of accuracy order two.
#' 3. Find the leftmost range of consecutive step sizes for which the slope of the truncation error
#'    is approximately equal (within 10% tolerance) to the accuracy order and which is long enough
#'    (e.g. at least length 5).
#' 4. Select the smallest step in this range and, if \code{correction = TRUE}, adjust
#'    it towards the predicted minimum of the combined rounding and truncation errors.
#'
#' The adjustment multiplies the selected step by
#' \eqn{r[(1-r^{-a})/(1+r^d)]^{1/(a+d)}}, where \eqn{r} is
#' \code{1/shrink.factor}, \eqn{d} is \code{deriv.order} and \eqn{a} is \code{acc.order}.
#'
#' Interpolation uses the existing derivative estimates and their leading dependence
#' on \code{h^acc.order}, avoiding further function evaluations. For continuous
#' functions, \eqn{f(x+\delta)\approx f(x)} for small \eqn{\delta}.
#' Only the existing grid is used; no higher-accuracy retries or extra candidate steps are evaluated.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical derivative (using central differences).
#'     \item \code{counts} – the number of paired grid points, with one additional central evaluation for even derivative orders.
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination due to a sequence of correct reductions.
#'         \item \code{1} – Reductions are slightly outside the tolerance.
#'         \item \code{2} – Tolerances are significantly violated; an approximate minimum is returned.
#'         \item \code{3} – Not enough finite function values; a rule-of-thumb value is returned.
#'         \item \code{4} – Step trimmed to 0.1|x| when |x| is not tiny and within range.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the step and argument grids,
#'       function values on those grids, estimated derivative values, and estimated error values.
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
#' f <- function(x) x^4  # The derivative at 1 is 4
#' step.M(x = 1, f)
#' step.M(x = 1, f, h0 = 1e-9) # Starting low
#' step.M(x = 1, f, h0 = 1000) # Starting high
#'
#' f <- sin  # The derivative at pi/4 is sqrt(2)/2
#' plot(step.M(x = pi/2, f))
#' plot(step.M(x = pi/4, f))
#' step.M(x = pi/4, f, h0 = 1e-9) # Starting low
#' step.M(x = pi/4, f, h0 = 1000) # Starting high
#' # where the truncation error estimate is invalid
step.M <- function(FUN, x, h0 = NULL, deriv.order = 1, acc.order = 2, range = NULL, shrink.factor = 0.5,
                   min.valid.slopes = 5L, seq.tol = 0.1, correction = TRUE,
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
  tstar <- (1 + inv.sf^deriv.order) / (1 - shrink.factor^acc.order)

  estimates <- differencesM(FUN, x, g, deriv.order, acc.order, inv.sf, retry = FALSE,
                            max.rel.error = max.rel.error, cores = cores, cl = cl,
                            preschedule = preschedule, ...)
  fgrid <- estimates$f
  cds <- estimates$deriv
  fhigher <- estimates$higher
  K <- estimates$K
  etrunc <- estimates$trunc
  eround <- estimates$round
  fallback <- NULL

  log2etrunc <- suppressWarnings(log2(etrunc))
  log2etrunc[!is.finite(log2etrunc)] <- NA
  dletrunc <- c(NA, diff(log2etrunc))
  dlh <- c(NA, diff(log2(hgrid)))
  slopes <- dletrunc / dlh

  # Prepaging an error message for the future
  sf.sugg <- max(0.5, round(sqrt(shrink.factor), 2))
  range.sugg <- range / c(1024, 1/1024)
  err1 <- paste0("Either increase 'shrink.factor' (e.g. from ", shrink.factor,
                 " to ", sf.sugg, ") to have a finer grid, or increase 'range' (",
                 "from [", toString(printE(range)), "] to ", "[", toString(printE(range.sugg)), "]).")

  # Find the first sequence of necessary length where the error is increasing
  hopt0 <- NA  # Overwritten in case of a successful check
  i.hopt <- i.increasing <- i.okay <- i.good <- NULL  # Overwritten in case of a successful check

  # Compare the observed slope with the target accuracy order
  good.slopes <- abs(slopes - acc.order) / acc.order <= seq.tol
  okay.slopes <- abs(slopes - acc.order) / acc.order <= min(seq.tol * 3, 0.5)
  positive.slopes <- slopes > 0

  # Now we need to run enough times before reaching the region dominated by rounding
  # If there is only one run, return its indices
  # Otherwise, return the indices of the first run of length 5 and more
  findFirstLongEnough <- function(x, min.len) {
    x[is.na(x)] <- FALSE
    if (sum(x, na.rm = TRUE) < 1) return(NULL)  # Nothing is finite = no runs
    r <- rle(x)
    tr <- which(r$values)
    if (max(r$lengths[tr]) < min.len) return(NULL)
    # Single run -- returning the indices without changes
    if (length(tr) == 1L) return(which(x))
    # Otherwise, find the first run of minimum length
    fr <- tr[which(r$lengths[tr] >= min.len)[1]]
    if (is.na(fr)) return(NULL)

    # Build the full-length logical vector with only that run set to TRUE
    out  <- rep(FALSE, length(x))
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    return(starts[fr]:ends[fr])
  }
  i.good <- findFirstLongEnough(good.slopes, min.len = min.valid.slopes)

  if (length(i.good) > 0) {
    i.hopt <- i.good[1]
    hopt0 <- hgrid[i.hopt]
  } else {  # Not enough good slopes -- use okay slopes
    i.okay <- findFirstLongEnough(okay.slopes, min.len = min.valid.slopes)
    if (length(i.okay) > 0) {
      i.hopt <- i.okay[1]
      hopt0 <- hgrid[i.hopt]
      exitcode <- 1L
      warning("The estimated truncation error has a slightly wrong reduction rate (~",
              round(stats::median(slopes[positive.slopes], na.rm = TRUE), 2), ", but should be ~", acc.order, "). ", err1)
    } else {
      if (any(is.finite(slopes))) {
        warning("The estimated truncation error has a wrong reduction rate (~", round(stats::median(slopes), 2),
                ", but should be ~", acc.order, "). ", err1)
      }
    }
  }

  i.increasing <- findFirstLongEnough(positive.slopes, min.len = max(length(i.good), length(i.okay), 2))

  if (is.finite(hopt0)) {
    # Finding the smallest step size with a valid slope and correcting it
    i.hopt <- which(hgrid == hopt0)
    # Mathur indexes the estimate at ratio*h, whereas our auxiliary stencil uses h
    hopt <- hopt0 * inv.sf * tstar^(-1/(deriv.order + acc.order))
  } else {  # Fail-safe return for cases like x^2 at 0
    fallback <- fallbackM(hgrid, fhigher, K, estimates$round.coef, eround,
                          x, deriv.order, acc.order)
    hopt0 <- hopt <- fallback$par
    i.hopt <- which.min(abs(log2(hgrid) - log2(hopt)))

    exitcode <- 2L
    if (sum(is.finite(fgrid)) < 3) exitcode <- 3L

    if (exitcode == 2L)
      warning(paste0("Could not find a sequence of ", min.valid.slopes, " reductions ",
                     "of the truncation error. ", err1,
                     " Finally, try setting 'min.valid.slopes' to 4 or even 3. For now, ",
                     "returning a very rough value based on the expected round-off error."))
    if (exitcode == 3L)
      warning(paste0("There are <3 finite function values on the grid. ",
                     "Try setting 'shrink.factor' to 0.9 (close to 1) or checking why the ",
                     "function does not return finite values on the range x+[",
                     toString(printE(range)), "] with ", n,
                     " exponentially spaced points. Returning a very rough value that may ",
                     "not even yield a finite numerical derivative."))
  }

  if (!correction) hopt <- hopt0
  # Final step size check
  if (hopt > 0.1*abs(x) && abs(x) > 4.71216091538e-7) {
    exitcode <- 4L
    h.target <- 0.1*abs(x)
    i.hopt <- which.min(abs(hgrid - h.target))
    hopt0 <- hopt <- hgrid[i.hopt]
  }

  msg <- switch(exitcode + 1L,
                "successfully found the optimal step size",  # 0
                "successfully found the optimal step size but allowed inaccurate slopes",  # 1
                "truncation branch ill-behaved, using a rounding-error fallback",  # 2
                "fewer than 3 finite function values on the grid, returning a rough step",  # 3
                "step size too large relative to x, using |x|/10 instead")  # 4

  diag.list <- list(h = hgrid, x = xgrid, f = fgrid, deriv = cds,
                    est.error = cbind(trunc = etrunc, round = eround),
                    slopes = slopes, max.rel.error = max.rel.error,
                    fhigher = fhigher, aprime = estimates$aprime,
                    fallback = fallback, f0 = estimates$f0)

  himin <- suppressWarnings(min(i.hopt, i.good, i.okay, i.increasing))
  i.round <- if (is.finite(himin)) 1:himin else NULL

  # Keep the undefined first slope aligned with the first grid point
  ret <- list(par = hopt, value = interpolateM(hgrid, cds, hopt, acc.order), counts = n,
              exitcode = exitcode, message = msg,
              abs.error = c(trunc = interpolateM(hgrid, etrunc, hopt, acc.order),
                            round = interpolateM(hgrid, eround, hopt, -deriv.order)),
              method = "Mathur", iterations = diag.list,
              args = list(deriv.order = deriv.order, acc.order = acc.order,
                          max.rel.error = max.rel.error, range = range, shrink.factor = shrink.factor,
                          min.valid.slopes = min.valid.slopes, seq.tol = seq.tol, correction = correction,
                          i.good = i.good, i.okay = i.okay, i.increasing = i.increasing, i.round = i.round))

  class(ret) <- "stepsize"
  return(ret)
}


plotM <- function(x, ...) {
  cols <- c("#7e1fde", "#328d2d", "#d58726", "#ca203a")
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
       main = "Mathur step-size selection", ...)
  # graphics::mtext(paste0("max. rel. err. = ", printE(x$args$max.rel.error, 1)), cex = 0.8, line = 0.5)
  if (length(x$args$i.round) > 0)
    graphics::points(h[x$args$i.round], et[x$args$i.round], pch = 16, col = cols[1], cex = 0.9)
  if (length(x$args$i.increasing) > 0)
    graphics::points(h[x$args$i.increasing], et[x$args$i.increasing], pch = 16, col = cols[4], cex = 0.9)
  if (length(x$args$i.okay) > 0)
    graphics::points(h[x$args$i.okay], et[x$args$i.okay], pch = 16, col = cols[3], cex = 0.9)
  if (length(x$args$i.good) > 0)
    graphics::points(h[x$args$i.good], et[x$args$i.good], pch = 16, col = cols[2], cex = 0.9)
  i.bad <- setdiff(seq_along(h), c(x$args$i.good, x$args$i.okay, x$args$i.increasing, x$args$i.round))
  if (length(i.bad) > 0)
    graphics::points(h[i.bad], et[i.bad], pch = 4, col = "#00000088", cex = 0.9)
  graphics::points(h[i.missing], rep(min(et, na.rm = TRUE), length(i.missing)),
                   pch = 4, col = "#00000088", cex = 0.9)

  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Rounding", "Trunc. good", "Trunc. fair", "Trunc. increas.", "Invalid"),
                   pch = c(16, 16, 16, 16, 4), col = c(cols, "#000000"),
                   box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(x)
}
