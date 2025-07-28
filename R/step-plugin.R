# Generator for plugin
getValsPlugin <- function(FUN, x, h, stage, cores, cl, preschedule, ...) {
  pow <- if (stage == 1) 3 else 1
  s <- fdCoef(deriv.order = pow)
  xgrid <- x + s$stencil * h * if (stage == 1) .Machine$double.eps^(-2/15) else 1

  FUNsafe <- function(z) safeF(FUN, z, ...)
  fp <- runParallel(FUN = FUNsafe, x = xgrid, cores = cores, cl = cl, preschedule = preschedule)
  has.errors <- any(e <- sapply(fp, function(y) inherits(attr(y, "error"), "error")))
  if (has.errors) {
    first.err.ind <- which(e)[1]
    warning(paste0("First error: ", as.character(attr(fp[[first.err.ind]], "error"))))
  }
  fgrid <- unlist(fp)

  f0 <- if (stage == 1) mean(fgrid[2:3]) else mean(fgrid)  # An approximation to f(x)
  cd <- sum(fgrid * s$weights) / h^pow

  list(x = xgrid, f = fgrid, cd = cd, f0 = f0)
}


#' Plug-in step selection
#'
#' @inheritParams step.DV
#'
#' @details
#' This function computes the optimal step size for central differences using the
#' plug-in approach.
#' The optimal step size is determined as the minimiser of the total error, which for central
#' finite differences is (assuming minimal bounds for relative rounding errors)
#' \deqn{\sqrt[3]{1.5 \frac{f'(x)}{f'''(x)} \epsilon_{\mathrm{mach}}}.}{(1.5 mach.eps * f' / f''')^(1/3).}
#' If the estimated third derivative is too small, the function assumes a third
#' derivative of 1 to prevent division-by-zero errors.
#'
#' @return A list similar to the one returned by \code{optim()}:
#'   \itemize{
#'     \item \code{par} – the optimal step size found.
#'     \item \code{value} – the estimated numerical first derivative (using central differences).
#'     \item \code{counts} – the number of iterations (in this case, it is 2).
#'     \item \code{abs.error} – an estimate of the truncation and rounding errors.
#'     \item \code{exitcode} – an integer code indicating the termination status:
#'       \itemize{
#'         \item \code{0} – Termination with checks passed tolerance.
#'         \item \code{1} – Third derivative is exactly zero; large step size preferred.
#'         \item \code{2} – Third derivative is too close to zero; large step size preferred.
#'         \item \code{3} – Solution lies at the boundary of the allowed value range.
#'         \item \code{4} – Step trimmed to 0.1|x| when |x| is not tiny and within range.
#'       }
#'     \item \code{message} – A summary message of the exit status.
#'     \item \code{iterations} – A list including the two-step size search path, argument grids,
#'       function values on those grids, and estimated third derivative values.
#'   }
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' f <- function(x) x^4
#' step.plugin(x = 2, f)
#' step.plugin(x = 0, f)  # f''' = 0, setting a large one
step.plugin <- function(FUN, x, h0 = max(1e-5*abs(x), stepx(x, deriv.order = 3)),
                        max.rel.error = .Machine$double.eps^(7/8), range = h0 / c(1e4, 1e-4),
                        cores = 1, preschedule = getOption("pnd.preschedule", TRUE),
                        cl = NULL, ...) {
  # TODO: add zero.tol everywhere
  if (length(x) != 1) stop(paste0("Direct step-size selection can handle only univariate inputs. ",
                                  "For 'x' longer than 1, use 'gradstep'."))
  cores <- checkCores(cores)
  h0 <- unname(h0)  # To prevent errors with derivative names
  cores <- min(cores, 4)
  if (length(range) != 2 || any(range <= 0)) stop("The range must be a positive vector of length 2.")
  range <- sort(range)

  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (inherits(cl, "cluster")) cores <- min(length(cl), cores)

  iters <- vector("list", 2)
  iters[[1]] <- getValsPlugin(FUN = FUN, x = x, h = h0, stage = 1,
                              cores = cores, cl = cl, preschedule = preschedule, ...)
  if (any(bad <- !is.finite(iters[[1]]$f))) {
    bad.iters <- 0
    while (TRUE) {
      bad.iters <- bad.iters + 1
      h0 <- h0 / 2
      if (h0 < max(range[1], .Machine$double.eps))
        stop(paste0("step.DV: Could not compute the function value at ", toString(iters[[1]]$x[bad]),
                    " after ", bad.iters, " attempts of step shrinkage",
                    ".\nChange the range, which is currently [", toString(range),
                    "], and/or\ntry a different starting h0, which is currently ", h0, "."))
      iters[[1]] <- getValsPlugin(FUN = FUN, x = x, h = h0, stage = 1,
                                  cores = cores, cl = cl, preschedule = preschedule, ...)
      if (!any(bad <- !is.finite(iters[[1]]$f))) break
    }
  }


  cd3 <- iters[[1]]$cd
  f0 <- iters[[1]]$f0

  exitcode <- 0L
  # If the estimate of f''' is near-zero, the step-size estimate may be too large --
  # only the modified one needs not be saved
  me13 <- max.rel.error^(1/3)
  if (abs(cd3) == 0) {
    exitcode <- 1L
    h <- pmax(me13, abs(x) / 128)
    cd3 <- me13^2 * abs(x)
  } else if (max(abs(f0), me13^2) / abs(cd3) > sqrt(1/max.rel.error)) {
    # The ratio of f' to f''' is too large -- safeguard against large steps
    # small values of f0 are truncated to macheps^(2/3) ~ 4e-11
    cd3 <- sqrt(max.rel.error) * max(abs(f0), me13^2)
    h <- pmax(me13, abs(x) / 256)
    exitcode <- 2L
  } else {  # Everything is OK
    h <- abs(1.5 * f0/cd3 * max.rel.error)^(1/3)
  }

  # x beyond the boundary
  if (h < range[1]) {
    h <- range[1]
    exitcode <- 3L
    side <- "left"
  } else if (h > range[2]) {
    h <- range[2]
    exitcode <- 3L
    side <- "right"
  }
  if (h > 0.1*abs(x) && abs(x) > 4.71216091538e-7) {
    # The found step size exceeds 10% of |x|, and x is not too small
    # FUN might poorly behave at x+h and x-h due to large steps.
    # Magic constant: threshold because stepx(4.712e-7) = 4.712e-8, i.e. 10%
    exitcode <- 4L
    h <- 0.1*abs(x)
  }

  iters[[2]] <- getValsPlugin(FUN = FUN, x = x, h = h, stage = 2,
                              cores = cores, cl = cl, preschedule = preschedule, ...)

  msg <- switch(exitcode + 1L,
                "successfully computed non-zero f''' and f'",  # 0
                "truncation error is zero, large step is favoured",  # 1
                "truncation error is near-zero, large step is favoured",  # 2
                paste0("step size too close to the ", side,
                       " end of the reasonable range [", printE(range[1]), ", ",
                       printE(range[2]), "]"),  # 3
                "step size too large relative to x, using |x|/10 instead")  # 4

  diag.list <- list(h = c(f3 = h0 * .Machine$double.eps^(-2/15), f1 = h),
                    x = unlist(lapply(iters, "[[", "x")),
                    f = unlist(lapply(iters, "[[", "f")),
                    deriv = c(cd3 = cd3, cd = iters[[2]]$cd),
                    args = list(h0 = h0, max.rel.error = max.rel.error, range = range))

  etrunc <- abs(cd3) / 6 * h^2
  eround <- 0.5 * .Machine$double.eps * max(abs(iters[[2]]$f)) / h
  ret <- list(par = h, value = iters[[2]]$cd, counts = 2, exitcode = exitcode,
              message = msg, abs.error = c(trunc = etrunc, round = eround),
              method = "plug-in", iterations = diag.list)
  class(ret) <- "stepsize"
  return(ret)

}

plotPlugin <- function(x, ...) {
  it <- x$iterations
  xl <- range(it$x)
  yl <- range(it$f)
  plot(it$x, it$f, bty = "n", xlim = xl, ylim = yl, col = c(rep("black", 4), rep("#FF0022", 2)),
       pch = 16, ylab = "Function values", xlab = "Step size", main = "Plug-in step-size selection", ...)
  graphics::abline(a = mean(it$f[5:6]) - it$deriv["cd"]*mean(it$x[5:6]), b = it$deriv["cd"], lty = 2)
  graphics::mtext(paste0("max. rel. err.: ", printE(it$args$max.rel.error, 1)), cex = 0.8, line = 0.5)
  return(invisible(x))
}

