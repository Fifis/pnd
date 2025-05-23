# Generator for Mathur and Kostyrka
getValsM <- function(FUN, x, cores, cl, preschedule, ...) {
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = x, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)
  return(matrix(fgrid, ncol = 2))
}


# Grid value generator for Mathur and Kostyrka
gridM <- function(x, range, shrink.factor = 0.5) {
  inv.sf <- 1/shrink.factor
  hgrid <- inv.sf^(floor(log(range[1], base = inv.sf)):ceiling(log(range[2], base = inv.sf)))
  xgrid   <- x + c(-hgrid, hgrid)
  list(h = hgrid, x = matrix(xgrid, ncol = 2L))
}


#' Mathur's AutoDX-like automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @inheritParams GenD
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param max.rel.error Error bound for the relative function-evaluation error
#'   (\eqn{\frac{\hat f(\hat x) - f(x)}{f(x)}}{(^f(^x) - f(x))/f(x)}). Measures how noisy a function is.
#'   If the function is relying on numerical optimisation routines, consider setting to
#'   \code{sqrt(.Machine$double.eps)}.
#'   If the function has full precision to the last bit, set to \code{.Machine$double.eps/2}.
#' @param shrink.factor A scalar less than 1 that is used to create a sequence of
#'   step sizes. The recommended value is 0.5. Change to 0.25 for a faster search. This
#'   number should be a negative power of 2 for the most accurate representation.
#' @param seq.tol Numeric scalar: maximum relative difference between old and new
#'   step sizes for declaring convergence.
#' @param min.valid.slopes Positive integer: how many points must form a sequence
#'   with the correct slope with relative difference from 2 less than \code{seq.tol}.
#'   If \code{shrink.factor} is small (< 0.33), consider reducing this to 4.
#' @param correction Logical: if \code{TRUE}, returns the corrected step size (last
#'   point in the sequence times a less-than-1 number to account for the possible
#'   continuation of the downwards slope of the total error); otherwise, returns
#'   the grid point that is is lowest in the increasing sequence of valid error
#'   estimates.
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{mathur2012analytical}{pnd} algorithm. It consists of the following steps.
#'
#' 1. Choose a reasonable large (but not too large) initial step size \eqn{h_0}{h0} and a reduction
#' factor (1/2 for fast, 1/4 for slow functions is a reasonable choice).
#' 2. Compute a series of truncation error estimates via third derivatives or Richardson extrapolation.
#' 3. Find the leftmost range of consecituve step sizes for which the slope of the trunctation error
#'    it approximately equal (within 10% tolerance) to the accuracy order and which is long enough
#'    (e.g. at least length 5).
#' 4. Use the leftmost point of this range as the uncorrected optimal step size, or correct it by shrinking
#'    it by a small amount given in the article.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (each iteration includes two function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination due to a sequence of correct reductions.
#'         \item \code{1} – Reductions are slightly outside the tolerance.
#'         \item \code{2} – Tolerances are significantly violated; an approximate minimum is returned.
#'         \item \code{3} – Not enough finite function values; a rule-of-thumb value is returned.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the step and argument grids,
#'       function values on those grids, estimated derivative values, and estimated error values.
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
  tstar <- (1 + inv.sf) / (1 - shrink.factor^2)  # TODO: n = 2...

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

  log2etrunc <- suppressWarnings(log2(etrunc))
  log2etrunc[!is.finite(log2etrunc)] <- NA
  dletrunc <- c(NA, diff(log2etrunc))
  dlh <- c(NA, diff(log2(hgrid)))
  slopes <- dletrunc / dlh

  # Rounding error for later
  stc <- fdCoef(deriv.order = deriv.order, acc.order = acc.order)
  l <- length(stc$stencil) / 2
  fmax <- apply(abs(cbind(fpos[, 1:l], fneg[, 1:l])), 1, max)
  f.eps <- max.rel.error * sum(abs(stc$weights)) * fmax
  f.delta <- 0.5 * .Machine$double.eps * fmax
  eround <- (f.eps + f.delta) / hgrid

  # Prepaging an error message for the future
  sf.sugg <- max(0.5, round(sqrt(shrink.factor), 2))
  range.sugg <- range / c(1024, 1/1024)
  err1 <- paste0("Either increase 'shrink.factor' (e.g. from ", shrink.factor,
                 " to ", sf.sugg, ") to have a finer grid, or increase 'range' (",
                 "from [", toString(printE(range)), "] to ", "[", toString(printE(range.sugg)), "]).")

  # Find the first sequence of necessary length where the error is increasing
  hopt0 <- NA  # Overwritten in case of a successful check
  i.hopt <- i.increasing <- i.okay <- i.good <- NULL  # Overwritten in case of a successful check

  # Checking the closeness to the slope to the approximation order, e.g. n=2
  good.slopes <- abs(slopes - acc.order) / acc.order <= seq.tol  # TODO: generalise with (d)
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
              round(stats::median(slopes[positive.slopes], na.rm = TRUE), 2), ", but should be ~2). ", err1)
    } else {
      if (any(is.finite(slopes))) {
        warning("The estimated truncation error has a wrong reduction rate (~", round(stats::median(slopes), 2),
                ", but should be ~2). ", err1)
      }
    }
  }

  i.increasing <- findFirstLongEnough(positive.slopes, min.len = max(length(i.good), length(i.okay), 2))

  if (is.finite(hopt0)) {
    # Finding the smallest step size with a valid slope and correcting it
    i.hopt <- which(hgrid == hopt0)
    hopt <- hopt0 * (1 / tstar)^(1/3)
  } else {  # Fail-safe return for cases like x^2 at 0
    # Two cases possible: f(x) = x^2 at x = 0 has eround increasing,
    # which implies that a large fail-safe step can be chosen
    f.er <- is.finite(eround)
    mean.sign <- if (sum(f.er) > 1) mean(sign(diff(eround[f.er]))) else 1
    # If the rounding error is growing, prefer a slightly larger step size
    if (mean.sign > 0.5) {
      hopt0 <- hopt <- 128*stepx(x = x, deriv.order = deriv.order, acc.order = acc.order,
                                 zero.tol = .Machine$double.eps^(1/3))
      i.hopt <- which.min(abs(log2(hgrid) - log2(hopt)))
    } else {
      # Otherwise, for an eround that does not grow convincingly,
      # find the step yielding the theoretical expected round-off value
      complete.rows <- apply(fgrid, 1, function(x) all(is.finite(x)))
      f0 <- abs(max(abs(fgrid[which(complete.rows)[1], ])))
      expected.eround <- (.Machine$double.eps^2 * f0^2 / 12)^(1/3)
      # TODO: generalise
      i.hopt <- which.min(abs(eround - expected.eround))
      hopt0 <- hopt <- hgrid[i.hopt]
    }

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

  msg <- switch(exitcode + 1L,
                "successfully found the optimal step size",
                "successfully found the optimal step size but allowed inaccurate slopes",
                "truncation error reduction rate is too wrong, returning the naive step",
                "Fewer than 3 finite function values on the grid, returning the naive step")

  diag.list <- list(h = hgrid, x = xgrid, f = fgrid, deriv = cds[[1]],
                    est.error = cbind(trunc = etrunc, round = eround),
                    slopes = slopes, max.rel.error = max.rel.error)

  himin <- suppressWarnings(min(i.hopt, i.good, i.okay, i.increasing))
  i.round <- if (is.finite(himin)) 1:himin else NULL

  # TODO: remove the first NA from the output
  ret <- list(par = if (correction) hopt else hopt0, value = cds[[1]][i.hopt], counts = n,
              exitcode = exitcode, message = msg, abs.error = diag.list$est.error[i.hopt, ],
              method = "Mathur", iterations = diag.list,
              args = list(max.rel.error = max.rel.error, range = range, shrink.factor = shrink.factor,
                          min.valid.slopes = min.valid.slopes, seq.tol = seq.tol, correction = correction,
                          i.good = i.good, i.okay = i.okay, i.increasing = i.increasing, i.round = i.round))

  class(ret) <- "stepsize"
  return(ret)
}


plotM <- function(x, ...) {
  cols <- c("#7e1fde", "#328d2d", "#d58726", "#ca203a")
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

  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Rounding", "Trunc. good", "Trunc. fair", "Trunc. increas.", "Invalid"),
                   pch = c(16, 16, 16, 16, 4), col = c(cols, "#000000"),
                   box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(x)
}
