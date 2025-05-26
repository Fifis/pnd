# Generator for Curtis--Reid
getValsCR <- function(FUN, x, h, deriv.order, acc.order, max.rel.error, cores, cl, preschedule, ...) {
  vanilla <- (deriv.order == 1) && (acc.order == 1)
  if (vanilla) {
    stc1 <- list(stencil = c(0, 1), weights = c(-1, 1))
    stc2 <- list(stencil = c(-1, 0, 1), weights = c(-0.5, 0, 0.5))
  } else {
    stc1  <- fdCoef(deriv.order = deriv.order, acc.order = acc.order)
    stc2  <- fdCoef(deriv.order = deriv.order, acc.order = acc.order+2)
  }

  xgrid <- x + stc2$stencil * h
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning(paste0("First error: ", as.character(attr(fp[[first.err.ind]], "error"))))
  }
  fgrid <- unlist(fp)
  n <- length(fgrid)

  if (vanilla) {
    less.acc <- sum(fgrid[-1] * stc1$weights) / h^deriv.order
  } else {
    idx <- match(stc1$stencil, stc2$stencil)
    if (!isTRUE(all.equal(idx, 2:(n-1))))
      stop("stepCR error: wrong indices, should never happen. Please file a bug report on GitHub.")
    less.acc <- sum(fgrid[idx] * stc1$weights) / h^deriv.order
  }
  more.acc <- sum(fgrid * stc2$weights) / h^deriv.order
  etrunc <- abs(less.acc - more.acc)
  eround <- max(abs(fgrid)) * max.rel.error / h^deriv.order
  if (identical(as.numeric(eround), 0)) eround <- .Machine$double.eps  # Safeguard against NA
  ratio <- etrunc / eround
  deriv <- if (vanilla) c(cd = more.acc, cd4 = NA) else c(cd = less.acc, cd4 = more.acc)
  ret <- list(h = h, x = xgrid, f = fgrid, ratio = ratio, deriv = deriv,
              est.error = c(trunc = etrunc, round = eround))
  return(ret)
}


#' Curtis--Reid automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @inheritParams fdCoef
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).
#' @param max.rel.error Error bound for the relative function-evaluation error
#'   (\eqn{\frac{\hat f(\hat x) - f(x)}{f(x)}}{(^f(^x) - f(x))/f(x)}). Measures how noisy a function is.
#'   If the function is relying on numerical optimisation routines, consider setting to
#'   \code{sqrt(.Machine$double.eps)}.
#'   If the function has full precision to the last bit, set to \code{.Machine$double.eps/2}.
#' @param aim Positive real scalar: desired ratio of truncation-to-rounding error. The original
#'   version over-estimates the truncation error, hence a higher \code{aim} is recommended.
#'   For the modernised version, aim should be equal to \code{deriv.order / acc.order}.
#' @param tol Numeric scalar greater than 1: tolerance multiplier for determining when to stop
#'   the algorithm based on the current estimate being between `aim/tol` and `aim*tol`.
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param maxit Integer: maximum number of algorithm iterations to prevent infinite
#'   loops in degenerate cases.
#' @param seq.tol Numeric scalar: maximum relative difference between old and new
#'   step sizes for declaring convergence.
#' @inheritParams runParallel
#' @param ... Passed to \code{FUN}.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{curtis1974choice}{pnd} algorithm.
#' If the estimated third derivative is exactly zero, then, the initial step size
#' is multiplied by 4 and returned.
#'
#' If 4th-order accuracy is requested, then, two things happen. Firstly,
#' since 4th-order-accurate differences requires a larger step size and the truncation error for
#' the 2nd-order-accurate differences grows if the step size is larger than the optimal one,
#' a higher ratio of truncation-to-rounding errors should be targeted. Secondly,
#' a 4th-order-accurate numerical derivative is returned, but the truncation and rounding errors
#' are still estimated for the 2nd-order-accurate differences. Therefore, the estimated truncation
#' error is higher and the real truncation error of 4OA differences is lower.
#'
#' TODO: mention that f must be one-dimensional
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences;
#'       especially useful for computationally expensive functions).
#'     \item \code{counts} – the number of iterations (each iteration includes three function evaluations).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Optimal termination within tolerance.
#'         \item \code{1} – Third derivative is zero; large step size preferred.
#'         \item \code{2} – No change in step size within tolerance.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Step trimmed to 0.1|x| when |x| is not tiny and within range.
#'         \item \code{5} – Maximum number of iterations reached.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the full step size search path, argument grids,
#'       function values on those grids, estimated error ratios, and estimated derivative values.
#'   }
#' @export
#'
#' @details
#' The arguments passed to \code{...} must not partially match those of [step.CR()]. For example, if
#' \code{cl} exists, then, attempting to avoid cluster export by using
#' \code{step.CR(f, x, h = 1e-4, cl = cl, a = a)} will result in an error: \code{a} matches \code{aim}
#' and \code{acc.order}. Redefine the function for this argument to have a name that is not equal
#' to the beginning of one of the arguments of [step.CR()].
#'
#' The original version of the article uses \code{deriv.order = 1}, \code{acc.order = 1}, \code{aim = 100}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.CR(x = 2, f)  # Wrap plot(...) around each of these functions
#' step.CR(x = 2, f, h0 = 1e-3)
#' step.CR(x = 2, f, acc.order = 2)
#' step.CR(x = 2, f, deriv.order = 2)
#' step.CR(x = 2, f, deriv.order = 3, acc.order = 4)
#'
#' # A bad start: too far away
#' step.CR(x = 2, f, h0 = 1000)  # Bad exit code + a suggestion to extend the range
#' step.CR(x = 2, f, h0 = 1000, range = c(1e-10, 1e5))  # Problem solved
#'
#' library(parallel)
#' cl <- makePSOCKcluster(names = 2, outfile = "")
#' abc <- 2
#' f <- function(x, abc) {Sys.sleep(0.02); abc*sin(x)}
#' x <- pi/4
#' system.time(step.CR(f, x, h = 1e-4, cores = 1, abc = abc))  # To remove speed-ups
#' system.time(step.CR(f, x, h = 1e-4, cores = 2, abc = abc))  # Faster
#' f2 <- function(x) f(x, abc)
#' clusterExport(cl, c("f2", "f", "abc"))
#' system.time(step.CR(f2, x, h = 1e-4, cl = cl))  # Also fast
#' stopCluster(cl)
step.CR <- function(FUN, x, deriv.order = 1, acc.order = 1, h0 = NULL,
                    max.rel.error = .Machine$double.eps^(7/8),
                    aim = NULL, tol = NULL, range = NULL, maxit = 20L, seq.tol = 1e-4,
                    cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                    cl = NULL, ...) {
  if (length(x) != 1) stop(paste0("Direct step-size selection can handle only univariate inputs. ",
                                  "For 'x' longer than 1, use 'gradstep'."))
  cores <- checkCores(cores)
  if (is.null(h0)) h0 <- stepx(x = x, deriv.order = deriv.order, acc.order = acc.order)
  h0 <- unname(h0)  # To prevent errors with derivative names
  deriv.order <- deriv.order[1]
  acc.order <- acc.order[1]
  vanilla <- (deriv.order == 1) && (acc.order == 1)
  if (deriv.order > 1 && acc.order == 1) acc.order <- 2
  cores <- min(cores, deriv.order+acc.order+1)
  if (is.null(range)) range <- h0 / c(1e5, 1e-5)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)
  if (range[1] < 2*.Machine$double.eps) range[1] <- 2*.Machine$double.eps  # Avoiding zeros
  if (is.null(tol)) tol <- if (vanilla) 10 else 4
  if (tol <= 1) stop("The tolerance 'tol' must be >=1 (e.g. 4).")
  if (is.null(aim)) aim <- if (vanilla) 100 else deriv.order / acc.order
  target <- sort(c(aim / tol, aim * tol))

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  iters <- list()
  i <- 1
  exitcode <- 0L

  while (i <= maxit) {
    hold <- if (i > 1) hnew else NA
    hnew <- if (i > 1) hold * (aim / max(iters[[i-1]]$ratio, 1))^(1/(deriv.order + acc.order)) else h0

    # Check the relative change of the step size, which is possible at
    # the 2nd iteration even before the 2nd error calculation
    if (i > 1) {
      if (abs(hnew/hold - 1) < seq.tol) {
        exitcode <- 2L
        break # Step 4: if the step size does not change, stop
      } else { # 4a, b: outside the range, replace with the border
        if (hnew < range[1]) hnew <- range[1]
        if (hnew > range[2]) hnew <- range[2]
        if (max(hnew/hold, hold/hnew) - 1 < seq.tol) {
          exitcode <- 3L
          break
        }
      }
    }

    res.i <- getValsCR(FUN = FUN, x = x, h = hnew, deriv.order = deriv.order,
                       acc.order = acc.order, max.rel.error = max.rel.error,
                       cores = cores, cl = cl, preschedule = preschedule, ...)
    iters[[i]] <- res.i
    if (any(bad <- !is.finite(res.i$f))) {
      bad.iters <- 0
      while (TRUE) {
        bad.iters <- bad.iters + 1
        hnew <- hnew / 2
        if (hnew < max(range[1], .Machine$double.eps))
          stop(paste0("step.CR: Could not compute the function value at ", toString(res.i$x[bad]),
                      " after ", bad.iters, " attempts of step shrinkage",
                      ".\nChange the range, which is currently [", toString(range),
                      "], and/or\ntry a different starting h0, which is currently ", h0, "."))
        res.i <- getValsCR(FUN = FUN, x = x, h = hnew, deriv.order = deriv.order,
                           acc.order = acc.order, max.rel.error = max.rel.error,
                           cores = cores, cl = cl, preschedule = preschedule, ...)
        if (!any(bad <- !is.finite(res.i$f))) break
      }
    }

    if (res.i$ratio >= target[1] && res.i$ratio <= target[2]) {
      break # Successful termination by matching the range
    }
    if (res.i$ratio < .Machine$double.eps^(4/5)) {
      exitcode <- 1L  # Zero truncation error -> only rounding error
      hnew <- hnew * 16 # For linear or quadratic functions, a large h is preferred
      i <- i+1
      iters[[i]] <- getValsCR(FUN = FUN, x = x, h = hnew, deriv.order = deriv.order,
                                acc.order = acc.order, max.rel.error = max.rel.error,
                                cores = cores, cl = cl, preschedule = preschedule, ...)
      break
    }

    i <- i + 1
  }

  i <- length(iters)
  if (i >= maxit) exitcode <- 5L

  if (hnew > 0.1*abs(x) && abs(x) > 4.71216091538e-7) {
    # The found step size exceeds 10% of |x|, and x is not too small
    # FUN might poorly behave at x+h and x-h due to large steps.
    # Magic constant: threshold because stepx(4.712e-7) = 4.712e-8, i.e. 10%
    exitcode <- 4L
    i <- i + 1
    hnew <- 0.1*abs(x)
    res.i <- getValsCR(FUN = FUN, x = x, h = hnew, deriv.order = deriv.order,
                              acc.order = acc.order, max.rel.error = max.rel.error,
                              cores = cores, cl = cl, preschedule = preschedule, ...)
    iters[[i]] <- res.i
  }

  msg <- switch(exitcode + 1L,
                "target error ratio reached within tolerance",  # 0
                "truncation error is close to zero, large step is favoured",  # 1
                "step size did not change between iterations",  # 2
                paste0("step size landed on the range ", if (res.i$h == range[1]) "left" else
                  "right", " end; consider extending the range"),  # 3
                "step size too large relative to x, using |x|/10 instead",  # 4
                "maximum number of iterations reached")  # 5

  diag.list <- list(h = do.call(c, lapply(iters, "[[", "h")),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    deriv = do.call(rbind, lapply(iters, "[[", "deriv")),
                    est.error = do.call(rbind, lapply(iters, "[[", "est.error")),
                    ratio = do.call(c, lapply(iters, "[[", "ratio")),
                    args = list(deriv.order = deriv.order, acc.order = acc.order,
                                h0 = h0, max.rel.error = max.rel.error, aim = aim,
                                tol = tol, range = range, maxit = maxit, seq.tol = seq.tol))

  ret <- list(par = res.i$h,
              value = unname(res.i$deriv["cd"]),
              counts = i, exitcode = exitcode, message = msg,
              abs.error = res.i$est.error,
              method = if (vanilla) "Curtis--Reid" else "Modified Curtis--Reid",
              iterations = diag.list)
  class(ret) <- "stepsize"
  return(ret)
}


plotCR <- function(x, ...) {
  it <- x$iterations
  et <- it$est.error[, "trunc"]
  er <- it$est.error[, "round"]
  et[et == 0] <- .Machine$double.eps
  er[er == 0] <- .Machine$double.eps
  cols <- c("#328d2d", "#7e1fde")
  xl <- it$args$range
  xl <- sqrt(xl * range(it$h))
  yl <- range(et[is.finite(et)], er[is.finite(er)]) * c(0.5, 2)
  vanilla <- (it$args$deriv.order == 1) && (it$args$acc.order == 1)
  plot(it$h, et, log = "xy", bty = "n", xlim = xl, ylim = yl,
       pch = 16, col = cols[1],
       ylab = "Estimated abs. error in df/dx", xlab = "Step size",
       main = paste0(if (!vanilla) "Modified " else NULL, "Curtis--Reid step-size selection"), ...)
  graphics::points(it$h, er, pch = 16, col = cols[2])
  if (length(it$h) > 1) {
    for (i in 2:length(it$h)) {
      graphics::arrows(it$h[i-1], et[i-1], it$h[i], et[i], angle = 20, length = 0.12, col = cols[1])
      graphics::arrows(it$h[i-1], er[i-1], it$h[i], er[i], angle = 20, length = 0.12, col = cols[2])
    }
  }
  ratios <- it$est.error[, "trunc"] / it$est.error[, "round"]
  ratios[ratios > 10] <- round(ratios[ratios > 10])
  ratios[1 < ratios & ratios < 10] <- round(ratios[1 < ratios & ratios < 10], 1)
  ratios[ratios < 1] <- round(ratios[ratios < 1], 3)
  graphics::text(it$h, et, labels = as.character(ratios), pos = ifelse(et > er, 1, 3))

  graphics::mtext(paste0("der. order ", it$args$deriv.order, ", acc. order ", it$args$acc.order,
                         ", target T/R ratio: [", round(it$args$aim/it$args$tol, 2),
                         ", ", round(it$args$aim*it$args$tol, 2),
                         "]; max. rel. err.: ", printE(it$args$max.rel.error, 1)),
                  cex = 0.8, line = 0.5)

  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Truncation", "Rounding"),
                   pch = 16, col = cols, box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 2)
  return(invisible(x))
}
