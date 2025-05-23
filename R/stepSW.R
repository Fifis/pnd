# Generator for Stepleman--Winarsky
getValsSW <- function(FUN, x, h, do.f0 = FALSE, ratio.last = NULL,
                      ratio.beforelast = NULL, max.rel.error,
                      cores, cl, preschedule, ...) {
  xgrid <- if (do.f0) x + c(-h, 0, h) else x + c(-h, h)

  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)

  cd <- (fgrid[length(fgrid)] - fgrid[1]) / (2*h)

  # Richardson extrapolation to estimate the truncation error
  etrunc <- if (is.null(ratio.last)) NA_real_ else abs((cd - ratio.last$deriv) / (1 - (ratio.last$h / h)^2))
  eround <- max.rel.error * max(abs(fgrid)) / h

  # Check monotonicity of the last three derivatives
  dlast <- if (!is.null(ratio.last) && !is.null(ratio.beforelast)) {
    c(ratio.beforelast$deriv, ratio.last$deriv, cd)
  } else {
    NULL
  }
  # Are derivative changes all positive or all negative?
  monotone.fp <- if (!is.null(dlast)) prod(diff(dlast)) > 0 else NA
  monotone.dfp <- if (!is.null(dlast)) abs(dlast[2] - dlast[3]) <= abs(dlast[2] - dlast[1]) else NA

  list(h = h, x = if (do.f0) xgrid[c(1, 3)] else xgrid, f = if (do.f0) fgrid[c(1, 3)] else fgrid,
       deriv = cd, f0 = if (do.f0) fgrid[2] else NULL, est.error = c(trunc = etrunc, round = eround),
       monotone = c(monotone.fp, monotone.dfp))
}


#' Stepleman--Winarsky automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param shrink.factor A scalar less than 1 that is used to multiply the step size
#'   during the search. The authors recommend 0.25, but this may be result in earlier
#'   termination at slightly sub-optimal steps. Change to 0.5 for a more thorough search.
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param seq.tol Numeric scalar: maximum relative difference between old and new
#'   step sizes for declaring convergence.
#' @param max.rel.error Positive numeric scalar > 0 indicating the maximum relative
#'   error of function evaluation. For highly accurate functions with all accurate bits
#'   is equal to half of machine epsilon. For noisy functions (derivatives, integrals,
#'   output of optimisation routines etc.), it is higher.
#' @param maxit Maximum number of algorithm iterations to avoid infinite loops.
#'   Consider trying some smaller or larger initial step size \code{h0}
#'   if this limit is reached.
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{stepleman1979adaptive}{pnd} algorithm.
#'
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (each iteration includes four function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination within tolerance.
#'         \item \code{2} – No change in step size within tolerance.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Maximum number of iterations reached.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the full step size search path, argument grids,
#'       function values on those grids, estimated derivative values, estimated error values,
#'       and monotonicity check results.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4  # The derivative at 1 is 4
#' step.SW(x = 1, f)
#' step.SW(x = 1, f, h0 = 1e-9) # Starting too low
#' # Starting somewhat high leads to too many preliminary iterations
#' step.SW(x = 1, f, h0 = 10)
#' step.SW(x = 1, f, h0 = 1000) # Starting absurdly high
#'
#' f <- sin  # The derivative at pi/4 is sqrt(2)/2
#' step.SW(x = pi/4, f)
#' step.SW(x = pi/4, f, h0 = 1e-9) # Starting too low
#' step.SW(x = pi/4, f, h0 = 0.1) # Starting slightly high
#' # The following two example fail because the truncation error estimate is invalid
#' step.SW(x = pi/4, f, h0 = 10)   # Warning
#' step.SW(x = pi/4, f, h0 = 1000) # Warning
step.SW <- function(FUN, x, h0 = 1e-5 * (abs(x) + (x == 0)),
                    shrink.factor = 0.5, range = h0 / c(1e12, 1e-8),
                    seq.tol = 1e-4, max.rel.error = .Machine$double.eps/2, maxit = 40L,
                    cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                    cl = NULL, ...) {
  if (length(x) != 1) stop(paste0("Direct step-size selection can handle only univariate inputs. ",
                                  "For 'x' longer than 1, use 'gradstep'."))
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 3)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  i <- 1
  exitcode <- 0L
  main.loop <- close.left <- FALSE
  first.main <- TRUE
  iters <- list()

  while (i <= maxit) {
    if (!main.loop) {
      if (i == 1) {
        res.i <- getValsSW(FUN = FUN, x = x, h = h0, max.rel.error = max.rel.error, do.f0 = TRUE,
                           ratio.last = NULL, ratio.beforelast = NULL,
                           cores = cores, cl = cl, preschedule = preschedule, ...)
        iters[[i]] <- res.i
        f0 <- res.i$f0
        hnew <- res.i$h
      }
      if (!is.finite(f0)) {
        stop(paste0("Could not compute the function value at ", x, ". FUN(x) must be finite."))
      }
      if (any(bad <- !is.finite(res.i$f))) {
        bad.iters <- 0
        while (TRUE) {
          bad.iters <- bad.iters + 1
          hnew <- hnew / 2
          if (hnew < max(range[1], .Machine$double.eps))
            stop(paste0("step.SW: Could not compute the function value at ", toString(res.i$x[bad]),
                        " after ", bad.iters, " attempts of step shrinkage",
                        ".\nChange the range, which is currently [", toString(range),
                        "], and/or\ntry a different starting h0, which is currently ", h0, "."))
          res.i <- getValsSW(FUN = FUN, x = x, h = h0, max.rel.error = max.rel.error, do.f0 = TRUE,
                             ratio.last = NULL, ratio.beforelast = NULL,
                             cores = cores, cl = cl, preschedule = preschedule, ...)
          if (!any(bad <- !is.finite(res.i$f))) break
        }
      }

      # First check: are the function values of different signs?
      # If yes, h is too large -- jump to the step-shrinking main loop:
      if (prod(sign(res.i$f)) < 0) { # f(x+h) should be close to f(x-h)
        main.loop <- TRUE
        i.prelim <- i
        next
      } else { # The bisection part
        while (!main.loop && i <= maxit) {
          hold <- hnew
          # N: approx. number of inaccurate digits due to rounding at h0 and hopt
          Nh0 <- log10(abs(f0)) - log10(2) - log10(hnew) - log10(abs(res.i$deriv))
          Nhopt <- -log10(max.rel.error) / 3
          rounding.nottoosmall <- Nh0 > 0 # Formula 3.16: some digits were lost,
          # the rounding error not zero, the step size not extremely large
          rounding.small <- Nh0 <= Nhopt + log10(shrink.factor) # Formula 3.15:
          # the rounding error is not worse than the truncation error
          if (rounding.small && rounding.nottoosmall) {
            main.loop <- TRUE # We are in the truncation branch -- go to shrinkage directly
            break
          } else { # Either rounding.small (too small) or rounding.nottoosmall (too large)
            i <- i + 1
            # This is what is necessary for 3.15 to hold with a small safety margin
            hnew <- hnew * 10^(Nh0 - (Nhopt + log10(shrink.factor)) + log10(2))
            # If 3.15 does not hold, this increases the step size
            # If 3.16 does not hold, shrink hnew towards the 3.15-optimal value
            if (!rounding.nottoosmall) hnew <- (hold + hnew) / 2 # Bisection

            if (hnew > range[2]) hnew <- range[2] # Safeguarding against huge steps
            res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                               ratio.last = iters[[i-1]], ratio.beforelast = NULL,
                               cores = cores, cl = cl, preschedule = preschedule, ...)
            iters[[i]] <- res.i
            if (abs(hnew/hold - 1) < seq.tol) { # The algorithm is stuck at one range end
              main.loop <- TRUE
              exitcode <- 2L
              break
            }

            bad <- !is.finite(res.i$f)
            if (any(bad) && !rounding.small) {  # Not in the original paper, but a necessary fail-safe
              warning(paste0("Could not compute the function value at [", toString(printE(res.i$x[bad])),
                             "]. FUN(", toString(printE(x)), ") is finite -- try the initial step h0 larger than ",
                             printE(h0), " but smaller than ", printE(hold), ". Halving from ",
                             printE(hnew), " to ", printE(hnew/2), ")."))
              for (i in 1:maxit) {
                hnew <- hnew/2
                res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                                   ratio.last = iters[[i-1]], ratio.beforelast = NULL,
                                   cores = cores, cl = cl, preschedule = preschedule, ...)
                iters[[i]] <- res.i
                if (is.finite(res.i$f)) break
              }
              if (!is.finite(res.i$f))
                stop(paste0("Could not compute the function value at [", toString(printE(res.i$x[bad])),
                            "]. FUN(", toString(printE(x)), ") is finite -- halving did not help.",
                            " Try 'gradstep(..., method = \"M\")' or 'step.M(...)' for a more reliable algorithm."))
            }

            if (any(bad) && !rounding.nottoosmall) {
              warning(paste0("Could not compute the function value at ", toString(res.i$x[bad, , drop = FALSE]),
                             ". FUN(", x, ") is finite -- try a step h0 smaller than ", hnew, ". ",
                             "Halving from ", printE(hnew), " to ", printE(hnew/2), ")."))
              for (i in 1:maxit) {
                hnew <- hnew/2
                res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                                   ratio.last = iters[[i-1]], ratio.beforelast = NULL,
                                   cores = cores, cl = cl, preschedule = preschedule, ...)
                iters[[i]] <- res.i
                if (is.finite(res.i$f)) break
              }
              if (!is.finite(res.i$f))
                stop(paste0("Could not compute the function value at [",
                            toString(printE(res.i$x[bad, , drop = FALSE])),
                            "]. FUN(", toString(printE(x)), ") is finite -- halving did not help.",
                            " Try 'gradstep(..., method = \"M\")' or 'step.M(...)' for a more reliable algorithm."))
            }
          }
        } # End initial step search
      }
      i.prelim <- i # Passing to the first iteration of the main loop with 2 function values saved
    } # End preliminary loop


    if (first.main) { # First main loop: extra h1 and h2 needed
      for (j in 1:2) { # Try a decreasing sequence
        i <- i + 1
        hnew <- hnew * shrink.factor
        if (hnew < range[1]) hnew <- range[1]
        res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                           ratio.last = iters[[i-1]], ratio.beforelast = if (j == 1) NULL else iters[[i-2]],
                           cores = cores, cl = cl, preschedule = preschedule, ...)

        iters[[i]] <- res.i
      }
      first.main <- FALSE
    } else { # Monotonicity satisfied, continuing shrinking, only h[i+1] needed
      if (abs(iters[[i]]$h/iters[[i-1]]$h - 1) < seq.tol) {
        exitcode <- 2L # If h did not shrink, it must have hit the lower bound
        break  # or something else went wrong; this code will most likely
        # be overwritten by 3; if it does not, throws a warning
      }

      hold <- hnew
      hnew <- hnew * shrink.factor
      if (hnew < range[1]) hnew <- range[1]
      i <- i + 1
      res.i <- getValsSW(FUN = FUN, x = x, h = hnew, max.rel.error = max.rel.error, do.f0 = FALSE,
                         ratio.last = iters[[i-1]], ratio.beforelast = iters[[i-2]],
                         cores = cores, cl = cl, preschedule = preschedule, ...)
      iters[[i]] <- res.i
    }

    if (any(!res.i$monotone)) break
  }
  hopt <- iters[[i-1]]$h # No monotonicity = bad
  hprev <- iters[[i-2]]$h
  # Error codes ordered by severity
  if (abs(hopt / hprev - 1) < seq.tol) exitcode <- 2L

  # If hprev hits the upper limit, the total error needs to be compared
  if (abs(hprev / range[2] - 1) < seq.tol) {
    abs.error.prev <- unname(iters[[i-1]]$est.error["trunc"] / shrink.factor + iters[[i-2]]$est.error["round"])
    abs.error.opt  <- unname(iters[[i-1]]$est.error["trunc"] + iters[[i-1]]$est.error["round"])
    if (abs.error.prev < abs.error.opt) { # The two-times reduction was unnecessary
      exitcode <- 3L
      hopt <- hprev
      iters[[i-2]]$est.error["trunc"] <- unname(iters[[i-1]]$est.error["trunc"] / shrink.factor)
    }
  }
  if (abs(hopt / range[1] - 1) < seq.tol) {
    exitcode <- 3L
    close.left <- TRUE
  }

  if (i >= maxit) exitcode <- 4L

  msg <- switch(exitcode + 1L,
                "successfully found a monotonicity violation",
                "", # The code cannot be 1 here
                "step size did not change between iterations",
                paste0("step size too close to the ", if (close.left)
                  "left" else "right", " end of the range; consider extending the range ",
                  "or starting from a ", if (close.left) "larger" else "smaller", " h0 value."),
                "maximum number of iterations reached")

  if (exitcode == 2L) warning(paste0("The step size did not change between iterations. ",
                                    "This should not happen. Send a bug report to https://github.com/Fifis/pnd/issues"))
  if (exitcode == 3L && !close.left)
    warning(paste0("The algorithm terminated at the right range of allowed step sizes. ",
                   "Possible reasons: (1) h0 is too low and the bisection step overshot ",
                   "the next value; (2) h0 was too large and the truncation error estimate ",
                   "is invalid; (3) the range is too narrow. Please try a slightly larger ",
                   "and a slightly smaller h0, or expand the range."))

  if (hopt > 0.1*abs(x) && abs(x) > sqrt(.Machine$double.eps)) {
    exitcode <- 3L
    warning(paste0("The found step size, ", printE(hopt), ", exceeds 10% of |x|, ",
                   abs(x), ", where x is not too small. FUN might poorly behave at x+h ",
                   "and x-h due to large steps. Try a different starting value h0 to be sure. ",
                   "Returning 0.01|x|."))
    i <- i + 1
    hopt <- 0.1*abs(x)
    res.i <- getValsSW(FUN = FUN, x = x, h = hopt, max.rel.error = max.rel.error,
                       do.f0 = FALSE, ratio.last = if (i > 1) iters[[i-1]] else NULL,
                       ratio.beforelast = if (i > 2) iters[[i-2]] else NULL,
                       cores = cores, cl = cl, preschedule = preschedule, ...)

    iters[[i]] <- res.i
  }

  diag.list <- list(h = do.call(c, lapply(iters, "[[", "h")),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    deriv = do.call(c, lapply(iters, "[[", "deriv")),
                    est.error = do.call(rbind, lapply(iters, "[[", "est.error")),
                    monotone = do.call(rbind, lapply(iters, "[[", "monotone")),
                    args = list(h0 = h0, shrink.factor = shrink.factor, range = range,
                                seq.tol = seq.tol, max.rel.error = max.rel.error, maxit = maxit))


  best.i <- if (exitcode == 3L && !close.left) i-2 else i-1
  ret <- list(par = hopt,
              value = iters[[best.i]]$deriv,
              counts = c(preliminary = i.prelim, main = i - i.prelim),
              exitcode = exitcode, message = msg,
              abs.error = iters[[best.i]]$est.error,
              method = "Stepleman--Winarsky",
              iterations = diag.list)
  class(ret) <- "stepsize"
  return(ret)
}


plotSW <- function(x, ...) {
  it <- x$iterations
  et <- it$est.error[, "trunc"]
  er <- it$est.error[, "round"]
  et[et == 0] <- .Machine$double.eps
  er[er == 0] <- .Machine$double.eps
  cols <- c("#328d2d", "#7e1fde")
  xl <- it$args$range
  xl <- sqrt(xl * range(it$h))
  yl <- range(et[is.finite(et)], er[is.finite(er)], na.rm = TRUE) * c(0.5, 2)
  plot(it$h[et > 0], et[et > 0], log = "xy", bty = "n", xlim = xl, ylim = yl,
       pch = 16, col = cols[1],
       ylab = "Estimated abs. error in df/dx", xlab = "Step size",
       main = paste0("Stepleman--Winarsky step-size selection"), ...)
  graphics::points(it$h, er, pch = 16, col = cols[2])
  if (length(it$h) > 1) {
    for (i in 2:length(it$h)) {
      graphics::arrows(it$h[i-1], et[i-1], it$h[i], et[i], angle = 20, length = 0.12, col = cols[1])
      graphics::arrows(it$h[i-1], er[i-1], it$h[i], er[i], angle = 20, length = 0.12, col = cols[2])
    }
  }

  graphics::mtext(paste0("max. rel. err.: ", printE(it$args$max.rel.error, 1)), cex = 0.8, line = 0.5)
  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Truncation", "Rounding"),
                   pch = 16, col = cols, box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(invisible(x))
}
