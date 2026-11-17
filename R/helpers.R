# Some of these internal functions are not exported
safeF <- function(FUN, x, ...) tryCatch(FUN(x, ...), error = function(e) return(structure(NA, error = e)))

# Keep derivative grids out of callbacks sent to workers; leave dots lazy.
bindSafeFunction <- function(FUN, ...) {
  force(FUN)
  function(y) safeF(FUN, y, ...)
}

# Resolve the function and arguments here, since socket workers have separate global environments
# A separate closure also avoids sending the caller's evaluation grid
bindFunctionArgs <- function(FUN, args) {
  force(FUN)
  force(args)
  function(x) do.call(FUN, c(list(x), args), quote = TRUE)
}

checkBadSafeF <- function(x) identical(as.logical(x), NA) && identical(names(attributes(x)), "error")

# Print in scientific (exponential) format like 1.23e-03 for 0.001234
printE <- function(x, d = 2) sprintf(paste0("%1.", d, ifelse(x >= 1/10^d & x <= 10^d, "f", "e")), x)

# Split a vector into contiguous runs
splitRuns <- function(x) {
  x <- sort(x)
  splits <- c(0, which(diff(x) > 1), length(x))
  runs <- vector("list", length(splits) - 1)
  for (i in seq_along(runs)) {
    start <- splits[i] + 1
    end <- splits[i + 1]
    runs[[i]] <- x[start:end]
  }
  runs
}

#' Print a matrix with separators
#'
#' @param x A numeric matrix to print line by line.
#' @inheritParams formatMat
#' @param begin A character to put at the beginning of each line, usually \code{""}, \code{"("}, or
#'   \code{"c("} (the latter is useful if console output is used in calculations).
#' @param sep The column delimiter, usually \code{"  "}, \code{"|"}, \code{"&"} (for LaTeX), or \code{", "}.
#' @param end A character to put at the end of each line, usually \code{""} or \code{")"}.
#' @param print If \code{TRUE}, outputs the lines of the matrix rows into the console.
#' @param format If \code{FALSE}, skips the formatting part.
#'
#' @returns An invisible character vector of formatted rows.
#' @export
#'
#' @examples
#' x <- matrix(c(-1234567, 12345.67, 123.4567,
#'               1.23456, -1.23456e-1, 0,
#'               1.23456e-4, 1.23456e-2, -1.23456e-6), nrow = 3)
#' printMat(x)
#' printMat(x, 2, TRUE, "c(", ", ", ")")  # Ready row vectors
printMat <- function(x, digits = 3, shave.spaces = TRUE,
                     begin = "", sep = "  ", end = "",
                     print = TRUE, format = TRUE) {
  if (format) x <- formatMat(x, digits = digits, shave.spaces = shave.spaces)
  x <- split(x, seq_len(nrow(x)))
  for (i in seq_along(x)) {
    x[[i]] <- paste0(begin, paste(x[[i]], collapse = sep), end, collapse = "")
    if (print) cat(x[[i]], "\n", sep = "")
  }
  return(invisible(unname(unlist(x))))
}

#' Round a matrix to N significant digits in mixed FP/exp notation
#'
#' @param x Numeric matrix.
#' @param digits Positive integer: the number of digits after the decimal point in scientific notation
#'   (i.e. one less than the number of significant digits).
#' @param shave.spaces Logical: if true, removes spaces to ensure compact output; if false, results
#'   in nearly fixed-width output (almost).
#'
#' @returns A character matrix of formatted entries.
#' @export
#'
#' @examples
#' x <- matrix(c(1234567, 12345.67, 123.4567,
#'               1.23456, -1.23456e-1, 0,
#'               -1.23456e-4, 1.23456e-2, -1.23456e-6), nrow = 3)
#' print(formatMat(x), quote = FALSE)
#' print(formatMat(x, digits = 1), quote = FALSE)
formatMat <- function(x, digits = 3, shave.spaces = TRUE) {
  econd <- abs(x) >= 0.1 & abs(x) < 10^(digits+1)
  fcond <- is.finite(x)
  xf <- x[econd & fcond]
  xe <- x[(!econd) & fcond]
  xi <- x[!fcond]

  if (length(xf) > 0) {  # Formatting like a float
    nd <- ceiling(log10(abs(xf)))  # Number of digits to the left of zero
    xfl <- split(xf, nd)
    if ("0" %in% names(xfl)) xfl[["0"]] <- sprintf(paste0("%1.", digits, "f"), xfl[["0"]])  # 0.12345
    for (i in 1:digits) {
      ii <- as.character(i)
      if (ii %in% names(xfl))
        xfl[[ii]] <- sprintf(paste0("%1.", digits+1-i, "f"), xfl[[ii]])
    }
    ii <- as.character(digits + 1)
    if (ii %in% names(xfl)) xfl[[ii]] <- paste0(as.character(round(xfl[[ii]])), ".")
    xf <- unsplit(xfl, nd)
  }

  if (length(xe) > 0) {  # Formatting in scientific format
    exact.zero <- which(xe == 0)
    xe <- sprintf(paste0("%1.", digits, "e"), xe)
    xe <- gsub("e([+-])0", "e\\1", xe)  # Shaving off the redundant zero
    if (any(exact.zero)) {
      xe[exact.zero] <- paste(c("0 ", rep(" ", digits)), collapse = "")
    }
  }

  xout <- vector("character", length(x))
  xout[econd & fcond] <- xf
  xout[(!econd) & fcond] <- xe

  has.minus <- grepl("^\\-", xout)  # Padding with spaces if there is a minus
  if (any(has.minus) && !all(has.minus)) xout[(!has.minus) & fcond] <- paste0(" ", xout[(!has.minus) & fcond])

  # Padding with spaces to right-align with the exponential caboose
  has.exp <- grepl("e", xout)
  if (any(has.exp) && !all(has.exp)) xout[(!has.exp) & fcond] <- paste0(xout[(!has.exp) & fcond],  "   ")

  # Restoring non-finite values
  xout[!fcond] <- format(xi, trim = TRUE)

  dim(xout) <- dim(x)
  if (is.null(dim(x))) xout <- matrix(xout, nrow = 1)

  # Removing extra spaces by column
  if (shave.spaces) {
    nhead <-  gsub("^( *).+", "\\1", xout)
    ntrail <- sub("^.*[^ ]", "", xout)
    nch <- nchar(xout)
    nh <- nchar(nhead)
    nt <- nchar(ntrail)
    hmin <- apply(nh, 2, min)
    tmin <- apply(nt, 2, min)
    # Removing redundancies by column
    for (i in which(hmin > 0)) {
      nchi <- nch[, i]
      nnew <- 1 + hmin[i]  # Being safe, not blindly substringing
      xout[, i] <- substr(xout[, i], nnew, nchi)
    }
    for (i in which(tmin > 0)) {
      nchi <- nch[, i]
      nnew <- nchi - tmin[i]  # Being safe, not blindly substringing
      xout[, i] <- substr(xout[, i], 1, nnew)
    }
  }

  return(xout)
}


#' Align printed output to the longest argument
#'
#' @param x A numeric vector or matrix to be aligned with a vector of column names.
#' @param names Optional: a character vector of element
#'   or column names to be output first. Overrides existing names.
#'   Numeric inputs are converted to character automatically.
#' @param pad A single character: \code{"l"} for left padding (flush-right justification),
#'   \code{"c"} for centre, and \code{"r"} for right padding (flush-left justification).
#'
#' @returns A character matrix with the first row of names and the rest aligned content
#' @export
#'
#' @examples
#' x <- structure(1:4, names = month.name[1:4])
#' print(alignStrings(x, names(x)), quote = FALSE)
#' print(alignStrings(x, names(x), pad = "c"), quote = FALSE)  # Centring
#' print(alignStrings(x, names(x), pad = "r"), quote = FALSE)  # Left alignment
#'
#' x <- matrix(c(1, 2.3, 4.567, 8, 9, 0), nrow = 2, byrow = TRUE)
#' colnames(x) <- c("Andy", "Bradley", "Ci")
#' alignStrings(x, pad = "c")
alignStrings <- function(x, names = NULL, pad = c("l", "c", "r")) {
  pad <- match.arg(pad)
  dimx <- dim(x)
  y <- if (!is.null(names)) as.character(names) else NULL
  if (is.null(y)) y <- if (is.null(dimx)) names(x) else colnames(x)
  if (is.null(dimx)) {
    x <- matrix(x, nrow = 1)
    dimx <- dim(x)
  }
  if (is.null(y)) {  # Nothing to align, exit
    x <- as.character(x)
    if (!is.null(dimx)) dim(x) <- dimx
    return(x)
  }
  if ((is.null(dimx) && length(x) != length(y)) || (!is.null(dimx) && ncol(x) != length(y)))
    stop("'x' must have the same length / number of columns as 'names'.")
  x <- as.character(x)
  if (!is.null(dimx)) dim(x) <- dimx
  yx <- rbind(y, x)
  n  <- nchar(yx)
  nmax <- apply(n, 2, max)
  nshort <- -sweep(n, 2, nmax, "-")
  repN <- function(n) if (n > 0) paste(rep(" ", n), collapse = "") else ""
  nleft <- switch(pad, l = nshort, c = floor(nshort/2), r = rep(0, length(nshort)))
  nright <- nshort - nleft
  pad.left  <- sapply(nleft, repN)
  pad.right <- sapply(nright, repN)
  yx <- paste0(pad.left, yx, pad.right)
  yx <- matrix(yx, nrow = if (is.null(dimx)) 2 else dimx[1]+1)
  return(yx)
}


#' Number of core checks and changes
#'
#' @param cores Integer specifying the number of CPU cores used for parallel computation.
#' Recommended to be set to the number of physical cores on the machine minus one.
#'
#' @returns An integer with the number of cores.
#' @export
#'
#' @examples
#' checkCores()
#' checkCores(2)
#' suppressWarnings(checkCores(1000))
checkCores <- function(cores = NULL) {
  if (identical(cores, 1) || identical(cores, 1L)) return(1)
  max.cores <- getOption("pnd.max.cores", 2L)
  if (is.null(cores)) cores <- max.cores
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != floor(cores)) stop("'cores' must be a positive integer.")

  if (cores > max.cores)
    warning("Requested more cores than you have logical cores. Consider setting 'cores = ", max.cores, "'.\n")

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")  # Limit to 2 cores for CRAN checks
  if (nzchar(chk) && chk == "TRUE") cores <- min(cores, 2L)
  return(cores)
}

#' Run a function in parallel over a list (internal use only)
#'
#' @param FUN A function of only one argument. If there are more arguments, use
#'   a wrapper such as \code{FUN2 <- function(x) FUN(x, a = 1)} and call it.
#' @param x A list to parallelise the evaluation of \code{FUN} over: either numbers
#'   or expressions.
#' @param preschedule Logical: if \code{TRUE}, enables pre-scheduling for \code{mclapply()}
#'   or disables load balancing with \code{parLapply()}. Recommended for functions that
#'   take less than 0.1 s per evaluation.
#' @param cores Integer specifying the number of CPU cores used for parallel computation.
#'   Recommended to be set to the number of physical cores on the machine minus one.
#' @param cl An optional user-supplied \code{cluster} object  (created by \code{makeCluster}
#'   or similar functions). If not \code{NULL}, the code uses \code{parLapply()} (if \code{preschedule}
#'   is \code{TRUE}) or \code{parLapplyLB()} on that cluster on all platforms.
#'
#'
#' @returns The value that `lapply(x, FUN)` would have returned.
#'
#' @details
#' If no custom cluster is supplied or no default cluster is registered, \code{options(pnd.parallel = "auto")}
#' uses \code{mclapply()} on Linux/macOS and automatically creates reusable PSOCK
#' workers on Windows. Set \code{options(pnd.parallel = "PSOCK")} to use these
#' workers on any platform, or \code{"mclapply"} to explicitly request forking on Linux/macOS.
#' Repeated calls reuse the pool; larger requests replace it and smaller requests
#' use a subset. Call [stopParallel()] to release it. Package unloading and session
#' exit also attempt clean-up. Supplied/default clusters remain owned by the caller;
#' their size overrides \code{cores} (the CRAN check limit of 2 is handled separately).
#' Nested calls inside automatically created workers run serially unless given a separate cluster.
#'
#' PSOCK workers are separate R sessions. The data needed by \code{FUN} should be in
#' its enclosing function environment or passed as arguments; objects used only
#' from the global environment need explicit export to a custom user-supplied cluster.
#' This also applies to unevaluated extra arguments referring to globals when
#' dimension checks are skipped. Use a local objective closure to avoid this dependency.
#' @export
#'
#' @examples
#' fslow <- function(x) Sys.sleep(x)
#' x <- rep(0.05, 6)
#' cl <- parallel::makeCluster(2)
#' print(t1 <- system.time(runParallel(fslow, x)))
#' print(t2 <- system.time(runParallel(fslow, x, cl = cl)))
#' print(t3 <- system.time(runParallel(fslow, x, cores = 2)))
#' parallel::stopCluster(cl)
#' cat("Parallel overhead at 2 cores: ", round(t2[3]*200/t1[3]-100), "%\n", sep = "")
#' # Ignore on Windows
#' cat("makeCluster() overhead at 2 cores: ", round(100*t2[3]/t3[3]-100), "%\n", sep = "")
#' # Store objective data locally so PSOCK workers can use it
#' f <- local({scale <- 2; function(x) scale * x^2})
#' runParallel(f, 1:3, cores = 2)
#' stopParallel()
runParallel <- function(FUN, x, cores = 1L, cl = NULL, preschedule = FALSE) {
  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (!is.null(cl)) {
    if (inherits(cl, "cluster"))
      cores <- length(cl)
    else  # Not a cluster implies rubbish input
      stop("The object passed as a cluster is not a cluster. ",
           "Use 'cl <- parallel::makeCluster(", cores, ")' to create a proper one.")
  }

  cores <- checkCores(cores)

  # If no specific parameters are requested, evaluate; if cl is not NULL, then, cores has the necessary length
  if (cores == 1 || !length(x) ||
      (is.null(cl) && (length(x) == 1L || isTRUE(getOption("pnd.parallel.worker", FALSE)))))
    return(lapply(x, FUN))

  if (is.null(cl)) cl <- getParallelCluster(cores)
  if (is.null(cl)) {
    return(parallel::mclapply(x, FUN, mc.cores = cores, mc.preschedule = preschedule))
  }

  if (inherits(cl, "cluster")) {
    cl <- cl[seq_len(min(cores, length(x)))]
    if (preschedule)
      return(parallel::parLapply(cl, x, FUN))
    else
      return(parallel::parLapplyLB(cl, x, FUN))
  }
  stop("'cl' should be 'lapply', 'mclapply', or a cluster object; 'cores' must be >= 1.")
}

#' Repeated indices of the first unique value
#'
#' @param m A matrix or a data frame.
#'
#' This function is an inverse function to such operations as
#' \code{m[c(1:3, 1, 1, 2), ]}: the matrix with potentially duplicated rows is
#' taken as input, and indices of each row in \code{unique(m)}
#' are returned.
#'
#' This function is faster -- at least in the examples tested so far -- than
#' `match(data.frame(t(m)), data.frame(t(unique(m))))`.
#'
#' @returns A vector of row indices in \code{unique(m)}.
#' @export
#'
#' @examples
#' dupRowInds(mtcars[rep(1:10, 10), rep(1:10, 10)])
#' dupRowInds(matrix(rnorm(1000), ncol = 10))
dupRowInds <- function(m) {
  if (is.data.frame(m)) m <- as.matrix(m)
  tm <- t(m)
  d <- duplicated(m)
  di <- which(d)
  if (length(di) == 0) return(seq_len(nrow(m)))
  mdup  <- m[di, , drop = FALSE]
  mudup  <- unique(mdup)  # Unique duplicates

  ret <- rep(NA, nrow(m))
  ret[!d] <- 1:sum(!d)
  for (i in seq_len(nrow(mudup))) {
    muniq.match <- which(apply(tm == mudup[i, ], 2, all))
    ret[muniq.match[-1]] <- ret[muniq.match[1]]
  }
  return(ret)
}


################################################################################
# The next block contains highly technical internal functions that exist only to
# outsource common routines from step.M() and step.K() into a somewhat modular
# code. Their outputs are most certainly not described in any of the
# accompanying articles, so use/change with caution. Both step.K and step.M
# attempt only one (preferably parallel) evaluation of the function on the
# grid because evaluating FUN is assumed to be more expensive that the
# overhead and the checks / algorithm themselves.

# The following function produces a b[j] = (0, +-1, +-2, +-4, ...)-like stencil.
# with weights w[j]. Exclude zero in case of symmetry (d+a-1 nodes suffice)
# Powers of ratio let wider stencils reuse function values from neighbouring trial steps
# offset records the shift: a node ratio^k*h_i uses the value at h_(i+k)
# Finally append the list element:
# $offset:  number k for which abs(b[j]) = r^k; NA for the centre b[j]=0.
stencilM <- function(d, a, ratio) {
  if (length(d) != 1L || !is.finite(d) || d < 1 || d != floor(d))
    stop("'deriv.order' must be a positive integer.")
  if (length(a) != 1L || !is.finite(a) || a < 2 || a %% 2L != 0)
    stop("'acc.order' must be a positive even integer.")
  positive <- ratio^(seq_len((d + a - 1L) %/% 2L) - 1L)
  near.int <- abs(positive/round(positive) - 1) < 16 * .Machine$double.eps
  positive[near.int] <- round(positive[near.int])
  nodes <- c(-rev(positive), if (d %% 2L == 0L) 0 else NULL, positive)
  stencil <- fdCoef(deriv.order = d, acc.order = a, stencil = nodes)
  stencil$offset <- match(abs(stencil$stencil), positive) - 1L
  stencil
}


# Construct the arguments at which FUN will be evaluated, once per grid point.
# The same geometric ratio must be used for candidate steps and stencil nodes
# This makes x +/- ratio^k*h_i available from the pair evaluated at h_(i+k)
# Columns of x hold the left and right arguments, with h increasing down the rows
# Start at the lower endpoint so a supplied range is not widened by rounding
gridM <- function(x, range, ratio = 2) {
  # Generate the next endpoint too, then discard it if it exceeds the upper bound
  # This avoids losing an admissible endpoint to rounding in the logarithms
  n <- 1L + ceiling((log(range[2]) - log(range[1])) / log(ratio))
  h <- range[1] * ratio^(seq_len(n) - 1L)
  # A wide interval can overflow ratio^i even when the resulting step is finite
  overflow <- !is.finite(h)
  h[overflow] <- exp(log(range[1]) + (which(overflow)-1L) * log(ratio))
  h <- h[is.finite(h) & h <= range[2]]
  list(h = h, x = matrix(x + c(-h, h), ncol = 2L))
}


# Function-value look-up: construct a table of previously calculated function values,
# arranging them in the order required by a finite-difference formula. One function
# evaluation can belong to several stencils.
# Each output row is a candidate step and each column corresponds to a node in stencil$stencil.
# Leave NA beyond the grid boundary rather than silently use an incomplete difference formula
stencilValuesM <- function(stencil, fgrid, f0, rows = seq_len(nrow(fgrid))) {
  values <- matrix(NA_real_, length(rows), length(stencil$weights))
  for (j in seq_along(stencil$weights)) {
    if (stencil$stencil[j] == 0) {
      values[, j] <- f0
    } else {
      # For the currently requested stencil position, find which row of the stored
      # function values contains the corresponding evaluation
      index <- rows + stencil$offset[j]
      # Larger steps may require points outside our grid
      inside <- index <= nrow(fgrid)
      # Negative stencil positions use the left column (1); ppsitive ones use the right column (2)
      values[inside, j] <- fgrid[cbind(index[inside], 1L + (stencil$stencil[j] > 0))]
    }
  }
  values
}

# Combine the assembled values into the weighted numerator of a difference.
# I.e. this function does NOT divide by h[i]^d; the caller does that.
# A central stencil is symmetric. For odd d, weights on the +/- positions
# have opposite signs, so use f(x+b*h)-f(x-b*h). For even d, weights have
# the same sign, so use their sum. Pairing first also allows equal values
# to cancel exactly for odd d, before multiplication by a weight.
stencilSumM <- function(values, stencil, d, check.zero = FALSE) {
  positive <- which(stencil$stencil > 0)
  negative <- match(-stencil$stencil[positive], stencil$stencil)
  if (d %% 2L == 1L) {
    pairs <- values[, positive, drop = FALSE] - values[, negative, drop = FALSE]
  } else {
    pairs <- values[, positive, drop = FALSE] + values[, negative, drop = FALSE]
  }
  result <- drop(pairs %*% stencil$weights[positive])
  centre <- which(stencil$stencil == 0)
  if (length(centre)) result <- result + values[, centre] * stencil$weights[centre]
  if (check.zero) {
    # A cancelled auxiliary sum can leave a tiny positive residue from fractional weights
    # Treat sums below their arithmetic resolution as unresolved, not as evidence for the fit
    # This checks the weighted sum of stored pairs, not the unknown error in FUN
    magnitude <- drop(abs(pairs) %*% abs(stencil$weights[positive]))
    if (length(centre)) magnitude <- magnitude + abs(values[, centre] * stencil$weights[centre])
    k <- 2L * length(positive) + length(centre)
    ku <- k * .Machine$double.eps/2
    unresolved <- is.finite(result) & abs(result) <= ku/(1-ku) * magnitude
    result[unresolved] <- 0
  }
  result
}

# Evaluate FUN once on a shared grid, then calculate derivatives and error
# model components for ALL candidate steps without further function calls.
# Estimating truncation error needs f^(d+a), whose stencil overlaps the target derivative's stencil
# A single batch of evaluations serves both; a higher-accuracy retry only recombines stored values.
# The key reuse idea: one stencil uses f(x +/- h[i]), another may use
# f(x +/- ratio*h[i]) = f(x +/- h[i+1]), which is already in the same cache.
differencesM <- function(FUN, x, grid, d, a, ratio, retry, max.rel.error,
                         cores, cl, preschedule, check.zero = FALSE, aprime = 2L,
                         richardson = FALSE, ...) {
  target <- stencilM(d, a, ratio)
  auxiliary <- stencilM(d + a, aprime, ratio)
  K <- abs(attr(target, "remainder.coef"))
  # Squared weights describe noise variance; shared nodes correlate adjacent derivatives
  # Use the primary auxiliary formula because it supplies the fitted rounding branch
  positive <- target$weights[target$stencil > 0]
  noise <- c(ratio = K * sqrt(sum(auxiliary$weights^2) / sum(target$weights^2)),
             overlap = (2 * sum(utils::head(positive, -1L) * utils::tail(positive, -1L)) +
                          sum(target$weights[target$stencil == 0]^2)) / sum(target$weights^2))
  n <- length(grid$h)
  centre <- d %% 2L == 0L
  f <- getValsM(FUN, c(grid$x, if (centre) x), cores, cl, preschedule, ...)
  fgrid <- matrix(f[seq_len(2L*n)], ncol = 2L)
  f0 <- if (centre) f[2L*n + 1L] else NA_real_
  values <- stencilValuesM(target, fgrid, f0)
  derivative <- stencilSumM(values, target, d) / grid$h^d
  fmax <- apply(abs(values), 1L, max)
  accuracy <- rep(aprime, n)
  if (richardson) {
    # Mathur subtracts the (already) rounded derivatives rather than regrouping their function values
    # The estimate belongs to the larger step of each adjacent pair
    # TODO: find page reference
    error <- c(NA_real_, diff(derivative)) / (1-ratio^(-a))
    higher <- error / attr(target, "remainder.coef")
    trunc <- abs(error)
  } else {
    values <- stencilValuesM(auxiliary, fgrid, f0)
    # d+a has the same parity as d because a is even, so the same pairing applies
    # Dividing by h^d leaves h^a*f^(d+a), the quantity needed for the truncation estimate
    higher <- stencilSumM(values, auxiliary, d, check.zero = check.zero) / grid$h^d
  }
  if (retry && !richardson) {
    # The shorter formula can still use inner values when an outer node is unavailable
    # Both formulae reuse the grid, so retrying adds no evaluations of FUN
    bad <- which(!is.finite(higher) | higher == 0)
    if (length(bad)) {
      alternative.accuracy <- if (aprime == 2L) 4L else 2L
      auxiliary <- stencilM(d + a, alternative.accuracy, ratio)
      values <- stencilValuesM(auxiliary, fgrid, f0, bad)
      alternative <- stencilSumM(values, auxiliary, d, check.zero = check.zero) / grid$h[bad]^d
      good <- is.finite(alternative) & alternative != 0
      higher[bad[good]] <- alternative[good]
      accuracy[bad[good]] <- alternative.accuracy
    }
  }
  if (!richardson) trunc <- K * abs(higher)
  round.coef <- (max.rel.error * sum(abs(target$weights)) + 0.5 * .Machine$double.eps) * fmax
  list(f = fgrid, f0 = f0, deriv = derivative, higher = higher / grid$h^a,
       trunc = trunc, round = round.coef / grid$h^d,
       K = K, aprime = accuracy, round.coef = round.coef, stencil = target$stencil, noise = noise)
}

# Choose the step for the interpolated derivative whose noise differs from a fresh stencil
# Assume independent errors of similar size in nearby function values
# The fitted branches give squared bias/variance; reused nodes supply their covariance
interpolationStepK <- function(h, usable, upper, gamma, d, a, noise) {
  lower <- utils::head(h, -1L)
  higher <- utils::tail(h, -1L)
  ratio <- higher/lower
  t <- ratio^a - 1
  variance <- ratio^(-2*d)
  covariance <- noise[["overlap"]] * ratio^(-d)
  # Divide squared bias by the variance at the lower step
  signal <- (noise[["ratio"]] * 2^((a+d)*(log2(lower)-gamma)))^2
  weight <- (1-covariance-signal*t) / (1+variance-2*covariance+signal*t^2)
  # Bound the interpolation weight before comparing pairs, including a partial pair at the cap
  maximum <- pmin(1, pmax(0, ((upper/lower)^a-1)/t))
  weight <- pmax(0, pmin(maximum, weight))
  step <- pmin(upper, lower * (1+t*weight)^(1/a))
  loss <- signal*(1+t*weight)^2 + (1-weight)^2 + variance*weight^2 +
    2*covariance*weight*(1-weight)
  # Different pairs have different variance scales; logarithms avoid overflow at distant steps
  loss <- log2(loss) - 2*d*(log2(lower)-gamma)
  complete <- utils::head(usable, -1L) & utils::tail(usable, -1L) & lower <= upper & is.finite(step) & step <= upper & is.finite(loss)
  if (!any(complete)) return(NA_real_)
  loss[!complete] <- Inf
  unname(step[which.min(loss)])
}

# For error A*h^a + B*h^(-d), the derivative of this expression gives an optimal rounding level
# Match that level on the observed grid instead of estimating a new derivative or calling FUN again
fallbackM <- function(h, higher, K, round.coef, eround, x, d, a) {
  finite.higher <- abs(higher[is.finite(higher)])
  # The median resists large noisy estimates at tiny steps; keep zeros rather than invent a positive coefficient
  A <- if (length(finite.higher)) K * stats::median(finite.higher) else NA_real_
  reference <- which(is.finite(round.coef) & round.coef > 0)
  # The finest complete stencil keeps the function-value scale local to x
  B <- if (length(reference)) round.coef[reference[1L]] else NA_real_
  # R* = B^(a/(a+d)) * (a*A/d)^(d/(a+d)), evaluated in logarithms
  target <- if (is.finite(A) && A > 0 && is.finite(B) && B > 0)
    exp((a * log(B) + d * (log(a) + log(A) - log(d))) / (a + d)) else NA_real_
  finite <- is.finite(eround)
  mean.sign <- if (sum(finite) > 1L) mean(sign(diff(eround[finite]))) else 1
  # Along increasing h, growth signals that the sampled function scale is not constant
  # Without ties, this threshold requires more than three quarters of differences to rise
  usable <- which(is.finite(eround) & eround > 0)
  if (mean.sign > 0.5 || !is.finite(target) || target <= 0 || !length(usable)) {
    # The factor 128 is an empirical escape from an unresolved constant-coefficient model
    par <- 128 * stepx(x, deriv.order = d, acc.order = a, zero.tol = .Machine$double.eps^(1/3))
  } else {
    par <- h[usable[which.min(abs(log(eround[usable]) - log(target)))]]
  }
  list(par = par, trunc.coef = A, round.coef = B, target = target, mean.sign = mean.sign)
}

# TODO: rewrite the following explanation in a less technical manner.
# When the truncation branch is ill-behaved, a larger step can reduce cancellation.
# Accept it only while the rounding intervals of all finer usable derivatives overlap.
# Three consecutive truncation estimates below their rounding bounds = choose the step
# using this rule. Also,
# compare all finer steps so an isolated coarse plateau is not accepted on its own.
stableStepM <- function(h, derivative, trunc, round, valid, upper, a) {
  if (is.complex(derivative)) return(integer(0))
  usable <- seq_along(h) %in% valid & is.finite(round)
  lower <- derivative - round
  higher <- derivative + round
  lower[!usable] <- -Inf
  higher[!usable] <- Inf
  compatible <- cummax(lower) <= cummin(higher)
  # A zero leading higher derivative can hide the next truncation term
  # Changes in target derivatives also estimate truncation, without extra function calls
  change <- c(abs(diff(derivative))/((h[-1L]/h[-length(h)])^a - 1), NA_real_)
  quiet <- usable & is.finite(trunc) & trunc >= 0 & trunc <= round & is.finite(change) & change <= round
  end <- seq_along(h)
  end <- end[end >= 3L]
  choices <- end[h[end] <= upper & compatible[end] & quiet[end-2L] & quiet[end-1L] & quiet[end]]
  if (length(choices)) choices[length(choices)] else integer(0)
}

# Avoid more function calls after step selection: the leading error varies linearly in h^a
# For the rounding bound use power=-d instead, preserving its expected inverse-power behaviour
interpolateM <- function(h, value, hopt, power) {
  closest <- sort(utils::head(order(abs(log(h) - log(hopt))), 2L))
  finite <- is.finite(value[closest])
  if (!any(finite)) {
    available <- which(is.finite(value))
    if (!length(available)) return(NA_real_)
    # Avoid the most extreme steps when neither neighbour is usable
    middle <- available[which.min(abs(available - length(h)/2))]
    return(unname(value[middle]))
  }
  closest <- closest[finite]
  if (length(closest) == 1L) return(unname(value[closest]))

  near <- h[closest]
  exact <- which(near == hopt)
  if (length(exact)) return(unname(value[closest[exact[1L]]]))
  # Reuse an endpoint instead of extrapolating beyond the two available steps
  if (hopt < min(near) || hopt > max(near))
    return(unname(value[closest[which.min(abs(log(near) - log(hopt)))]]))

  # Normalise so the powered coordinates do not exceed one, avoiding overflow
  scale <- if (power > 0) max(near) else min(near)
  nodes <- (near/scale)^power
  weight <- ((hopt/scale)^power - nodes[1L]) / diff(nodes)
  # Preserve equal derivatives exactly, rather than round two weighted copies separately
  difference <- value[closest[2L]] - value[closest[1L]]
  if (is.finite(difference)) return(unname(value[closest[1L]] + weight * difference))
  # Opposite large values can overflow their difference even when the weighted sum is finite
  unname(sum(value[closest] * c(1 - weight, weight)))
}


