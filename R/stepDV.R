# Generator for Dumontet--Vignes
getValsDV <- function(FUN, x, k, max.rel.error, cores, cl, preschedule, ...) {
  xgrid <- x + c(-2*k, -k, k, 2*k)
  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning("First error: ", as.character(attr(fp[[first.err.ind]], "error")))
  }
  fgrid <- unlist(fp)

  tgrid <- fgrid * c(-0.5, 1, -1, 0.5) # T1, ..., T4 from the paper
  A <- sum(tgrid[tgrid > 0]) # Only positive terms
  B <- sum(tgrid[tgrid < 0]) # Only negative terms
  # TODO: error if fgrid has different sign
  f3sup <- (A / (1-max.rel.error) + B / (1+max.rel.error)) / k^3
  f3inf <- (A / (1+max.rel.error) + B / (1-max.rel.error)) / k^3
  f3 <- sum(tgrid) / k^3 # Estimate of third derivative
  L <- f3sup / f3inf
  L <- max(abs(L), abs(1/L))
  ret <- list(k = k, x = xgrid, f = fgrid, ratio = L,
              deriv = c(f3inf = f3inf, f3 = f3, f3sup = f3sup))
  return(ret)
}


#' Dumontet--Vignes automatic step selection
#'
#' @param x Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.
#' @param FUN Function for which the optimal numerical derivative step size is needed.
#' @param h0 Numeric scalar: initial step size, defaulting to a relative step of
#'   slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}). This step
#'   size for first derivarives is internallt translated into the initial step size for third
#'   derivatives by multiplying it by the machine epsilon raised to the power -2/15.
#' @param range Numeric vector of length 2 defining the valid search range for the step size.
#' @param max.rel.error Positive numeric scalar > 0 indicating the maximum relative
#'   error of function evaluation. For highly accurate functions with all accurate bits
#'   is equal to \code{.Machine$double.eps/2}. For noisy functions (derivatives, integrals,
#'   output of optimisation routines etc.), it is higher, typically
#'   \code{sqrt(.Machine$double.eps)}. Dumontet and Vignes recommend
#'   \code{.Machine$double.eps^(3/4) = 2e-12} for common functions.
#' @param ratio.limits Numeric vector of length 2 defining the acceptable ranges
#'   for step size: the algorithm stops if the relative perturbation of the third derivative by
#'   amplified rounding errors falls within this range.
#' @param maxit Maximum number of algorithm iterations to avoid infinite loops in cases
#'   the desired relative perturbation factor cannot be achieved within the given \code{range}.
#'   Consider extending the range if this limit is reached.
#' @inheritParams runParallel
#' @param ... Passed to FUN.
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' \insertCite{dumontet1977determination}{pnd} algorithm.
#' If the estimated third derivative is exactly zero, the function assumes a third
#' derivative of 1 to prevent division-by-zero errors.
#'
#' Note: the iteration history tracks the third derivative, not the first.
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
#'         \item \code{1} – Third derivative is zero; large step size preferred.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Step trimmed to 0.1|x| when |x| is not tiny and within range.
#'         \item \code{5} – Maximum number of iterations reached; optimal step size is within the allowed range.
#'         \item \code{6} – Maximum number of iterations reached; optimal step size
#'           was outside allowed range and had to be snapped to a boundary or to 0.1|x|.
#'         \item \code{7} – No search was performed (used when \code{maxit = 1}).
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the full step size search path (note: for the third derivative),
#'       argument grids, function values on those grids, and estimated third derivative values.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.DV(x = 2, f)
#' step.DV(x = 2, f, h0 = 1e-3)
#'
#' # Alternative plug-in estimator with only one evaluation of f'''
#' step.DV(x = 2, f, maxit = 1)
#' step.plugin(x = 2, f)
step.DV <- function(FUN, x, h0 = 1e-5*max(abs(x), sqrt(.Machine$double.eps)),
                    range = h0 / c(1e6, 1e-6), max.rel.error = .Machine$double.eps^(7/8),
                    ratio.limits = c(2, 15), maxit = 40L,
                    cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                    cl = NULL, ...) {
  if (length(x) != 1) stop(paste0("Direct step-size selection can handle only univariate inputs. ",
                                  "For 'x' longer than 1, use 'gradstep'."))
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 4)
  k0 <- h0 * .Machine$double.eps^(-2/15)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)
  range3 <- range * .Machine$double.eps^(-2/15)

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  iters <- list()
  i <- 1
  exitcode <- 0L
  ndownwards <- 0L # For counting the number of downwards shrinkages

  while (i <= maxit) {
    if (i == 1) k <- k0
    res.i <- getValsDV(FUN = FUN, x = x, k = k, max.rel.error = max.rel.error,
                       cores = cores, cl = cl, preschedule = preschedule, ...)
    iters[[i]] <- res.i
    if (any(bad <- !is.finite(res.i$f))) {
      bad.iters <- 0
      while (TRUE) {
        bad.iters <- bad.iters + 1
        k <- k / 2
        if (k < max(range[1], .Machine$double.eps))
          stop(paste0("step.DV: Could not compute the function value at ", toString(res.i$x[bad]),
                      " after ", bad.iters, " attempts of step shrinkage",
                      ".\nChange the range, which is currently [", toString(range),
                      "], and/or\ntry a different starting h0, which is currently ", h0, "."))
        res.i <- suppressWarnings(getValsDV(FUN = FUN, x = x, k = k, max.rel.error = max.rel.error,
                                            cores = cores, cl = cl, preschedule = preschedule, ...))
        if (!any(bad <- !is.finite(res.i$f))) break
      }
    }

    # Quick rule of thumb: stop after the first iteration
    if (maxit == 1) {
      exitcode <- 7L
      break
    }

    # If the estimate of f''' is near-zero, then, the algorithm must stop to avoid division by near-0
    if (abs(res.i$deriv["f3"]) < 8 * .Machine$double.eps) {
      exitcode <- 1L
      break
    }

    # If the function sign fluctuates, the step size is too wide
    f.opposite.sign    <- !all(diff(sign(res.i$f)) == c(0, 0, 0))
    # If the numerically contaminated derivative estimates blow up, the step size is too small
    flim.opposite.sign <- !all(diff(sign(res.i$deriv)) == c(0, 0))

    if (res.i$ratio > ratio.limits[2] || flim.opposite.sign) {
      # The rounding error is too high; the step size must be increased
      range3[1] <- k
    } else {
      if (res.i$ratio < ratio.limits[1] || f.opposite.sign) {
        # The rounding error is too small; the step size must be decreased
        range3[2] <- k
        ndownwards <- ndownwards + 1
      } else {
        # The numerical error is in the good range, the step size is optimal
        break
      }
    }
    k <- sqrt(prod(range3)) # The interval has been sub-divided; take its geometric mean
    i <- i+1
  }
  i <- length(iters)

  f3 <- if (exitcode != 1L) sum(res.i$f * c(-0.5, 1, -1, 0.5)) / k^3 else 1
  f0 <- mean(res.i$f[2:3]) # Approximately f(x); the error is small for small h
  h <- (1.68 * max.rel.error * abs(f0/f3))^(1/3) # Formula 36 from Dumontet & Vignes (1977)

  if (h < range[1]) {
    h <- range[1]
    exitcode <- 3L
    side <- "left"
  }
  if (h > range[2]) {
    h <- range[2]
    exitcode <- 3L
    side <- "right"
  }

  if (h > 0.1*abs(x) && abs(x) > 4.71216091538e-7) {
    # The found step size exceeds 10% of |x|, and x is not too small
    # FUN might poorly behave at x+h and x-h due to large steps.
    # Magic constant: threshold because stepx(4.712e-7) = 4.712e-8, i.e. 10%
    # No iteration update because h is merely a by-product of optimisation
    h <- 0.1*abs(x)
    exitcode <- 4L
  }

  # Was the procedure systematically unsuccsessful?
  # This error code is more severe than
  if (i >= maxit && maxit > 1) {  # Did it waste many iterations in vain?
    exitcode <- if (h <= range[1] || h >= range[2]) 6L else 5L
    side <- if (ndownwards >= maxit/2) "right" else "left"
  }

  h <- h + x # Minimising the representation error
  h <- h - x
  xgrid <- c(x-h, x+h)
  fgrid <- vapply(xgrid, FUN, ..., FUN.VALUE = numeric(1))
  cd <- (fgrid[2] - fgrid[1]) / h / 2
  etrunc <- unname(abs(res.i$deriv["f3"])) * h^2 / 6     # Formula for 'em' from Dumontet (1973), 1.4.2 (p. 37)
  eround <- max.rel.error * max(abs(fgrid)) / h  # Formula for 'ed' ibid.

  msg <- switch(exitcode + 1L,
                "target error ratio reached within tolerance",  # 0
                "truncation error is zero, large step is favoured",  # 1
                "",  # 2 -- not used yet
                paste0("step size too close to the ", side,
                       " end of the range [", printE(range[1]), ", ",
                       printE(range[2]), "]; consider extending it"),  # 3
                "maximum number of iterations reached",  # 4
                paste0("maximum number of iterations reached and step size occured on the ",
                       side, " end of the range [", printE(range[1]), ", ",
                       printE(range[2]), "]; consider expanding it"),  # 5
                "only one iteration requested; rough values returned")  # 6

  diag.list <- list(k = do.call(c, lapply(iters, "[[", "k")),
                    x = do.call(rbind, lapply(iters, "[[", "x")),
                    f = do.call(rbind, lapply(iters, "[[", "f")),
                    ratio = do.call(c, lapply(iters, "[[", "ratio")),
                    deriv3 = do.call(rbind, lapply(iters, "[[", "deriv")),
                    args = list(h0 = h0, range = range, max.rel.error = max.rel.error,
                                ratio.limits = ratio.limits, maxit = maxit))

  ret <- list(par = h, value = cd, counts = i, exitcode = exitcode,
              message = msg, abs.error = c(trunc = etrunc, round = eround),
              method = "Dumontet--Vignes",
              iterations = diag.list)
  class(ret) <- "stepsize"
  return(ret)
}


plotDV <- function(x, ...) {
  it <- x$iterations
  f3inf <- it$deriv3[, "f3inf"]
  f3sup <- it$deriv3[, "f3sup"]
  f3 <- it$deriv3[, "f3"]
  f3inf <- sign(f3inf) * log1p(abs(f3inf))
  f3sup <- sign(f3sup) * log1p(abs(f3sup))
  f3 <- sign(f3) * log1p(abs(f3))
  cols <- c("#EE0022", "#4422BB", "#009900AA")
  xl <- sqrt(it$args$range * range(it$k))
  yl <- range(f3inf, f3sup)
  plot(it$k, f3inf, log = "x", bty = "n", xlim = xl, ylim = yl, yaxt = "n",
       pch = ifelse(f3inf > 0, 16, 0), col = cols[1],
       ylab = "Est. f''' and its bounds", xlab = "Step size",
       main = "Dumontet--Vignes step-size selection", ...)
  graphics::points(it$k, f3sup, pch = ifelse(f3sup > 0, 16, 0), col = cols[2])
  graphics::abline(h = 0, v = x$par, lty = 2)

  graphics::points(it$k, f3, pch = 1, cex = 0.8, col = cols[3])
  if (length(it$k) > 1) {
    for (i in 2:length(it$k)) graphics::arrows(it$k[i-1], f3[i-1], it$k[i], f3[i], angle = 20, length = 0.12, col = cols[3])
  }
  labs <- pretty(yl)
  labs <- sign(labs) * expm1(abs(labs))
  ratios <- sprintf("%1.1f", it$ratio)
  ratios[sign(f3inf) * sign(f3sup) < 0] <- "X"
  graphics::text(it$k, f3sup, labels = ratios, pos = c(1, 3))
  graphics::axis(2, pretty(yl), sprintf("%1.1e", labs))

  graphics::mtext(paste0("target sup/inf ratio: [", toString(round(it$args$ratio.limits, 1)),
                         "]; max. rel. err.: ", printE(it$args$max.rel.error, 1)),
                  cex = 0.8, line = 0.5)

  graphics::abline(v = x$par, lty = 3, col = "#00000088")
  graphics::legend("topleft", c("Lower", "Upper", "f'''"),
                   pch = 16, col = cols, box.col = "#FFFFFF00", bg = "#FFFFFFAA", ncol = 1)
  return(invisible(x))
}
