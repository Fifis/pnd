# Some of these internal functions are not exported

safeF <- function(FUN, x, ...) tryCatch(FUN(x, ...), error = function(e) return(structure(NA, error = e)))

checkBadSafeF <- function(x) identical(as.logical(x), NA) && identical(names(attributes(x)), "error")

# Concatenate together with a comma between the terms
pasteAnd <- function(x) paste0(x, collapse = ", ")

# Print in scientific (exponential) format like 1.23e-03 for 0.001234
printE <- function(x, d = 2) sprintf(paste0("%1.", d, ifelse(x >= 0.01 & x <= 10^d, "f", "e")), x)

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
#' @returns The same \code{x} that was passed as the first input.
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
    x[[i]] <- paste0(begin, paste0(x[[i]], collapse = sep), end, collapse = "")
    if (print) cat(x[[i]], "\n", sep = "")
  }
  return(invisible(unname(unlist(x))))
}

#' Round a matrix to N signifcant digits in mixed FP/exp notation
#'
#' @param x Numeric matrix.
#' @param digits Positive integer: the number of digits after the decimal comma to round to
#'   (i.e. one less than the number of significant digits).
#' @param shave.spaces Logical: if true, removes spaces to ensure compact output; if false, results
#'   in nearly fixed-width output (almost).
#'
#' @returns A numeric matrix with all entries of equal width with the same number of characters
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
  xi <- x[(!econd) & (!fcond)]

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
      xe[exact.zero] <- paste0(c("0 ", rep(" ", digits)), collapse = "")
    }
  }

  xout <- vector("character", length(x))
  xout[econd & fcond] <- xf
  xout[(!econd) & fcond] <- xe

  has.minus <- grepl("^\\-", xout)  # Padding with spaces if there is a minus
  if (any(has.minus) && any(!has.minus)) xout[(!has.minus) & fcond] <- paste0(" ", xout[(!has.minus) & fcond])

  # Padding with spaces to right-align with the exponential caboose
  has.exp <- grepl("e", xout)
  if (any(has.exp) && any(!has.exp)) xout[(!has.exp) & fcond] <- paste0(xout[(!has.exp) & fcond],  "   ")

  # Restoring non-finite values
  x[(!econd) & (!fcond)] <- xi

  dim(xout) <- dim(x)
  if (is.null(dim(x))) xout <- matrix(xout, nrow = 1)

  # Removing extra spaces by column
  if (shave.spaces) {
    nhead <-  gsub("^( *).+", "\\1", xout)
    ntrail <- gsub("^ *[.0-9e+-]+", "", xout)
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
#' @param names Optional: if x does not have (column) names, a character vector of element
#'   or column names to be output first. Ignored if \code{x} is named.
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
  repN <- function(n) if (n > 0) paste0(rep(" ", n), collapse = "") else ""
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

  if (cores > max.cores)
    warning(paste0("Requested more cores than you have logical cores. Consider setting 'cores = ", max.cores, "'.\n"))

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")  # Limit to 2 cores for CRAN checks
  if (nzchar(chk) && chk == "TRUE") cores <- min(cores, 2L)
  return(cores)
}

#' Run a function in parallel over a list (internal use only)
#'
#' @param FUN A function of only one argument. If there are more arguments, use
#'   the \code{FUN2 <- do.call(FUN, c(list(x), ...))} annd call it.
#' @param x A list to parallelise the evaluation of \code{FUN} over: either numbers
#'   or expressions.
#' @param preschedule Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
#'   or enables load balancing with \code{parLapplyLB()}. Recommended for functions that
#'   take less than 0.1 s per evaluation.
#' @param cores Integer specifying the number of CPU cores used for parallel computation.
#'   Recommended to be set to the number of physical cores on the machine minus one.
#' @param cl An optional user-supplied \code{cluster} object  (created by \code{makeCluster}
#'   or similar functions). If not \code{NULL}, the code uses \code{parLapply()} (if \code{preschedule}
#'   is \code{TRUE}) or \code{parLapplyLB()} on that cluster on Windows, and \code{mclapply}
#'   (fork cluster) on everything else.
#'
#'
#' @returns The value that `lapply(x, FUN)` would have returned.
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
runParallel <- function(FUN, x, cores = 1L, cl = NULL, preschedule = FALSE) {

  if (!is.null(cl)) {
    if (inherits(cl, "cluster"))
      cores <- length(cl)
    else  # Not a cluster implies rubbish input
      stop(paste0("The object passed as a cluster is not a cluster. ",
                  "Use 'cl <- parallel::makeCluster(", cores, ")' to create a proper one."))
  }

  cores <- checkCores(cores)

  # If no specific parameters are requested, evaluate; if cl is not NULL, then, cores has the necessary length
  if (cores == 1) {
    #  if (is.list(x) && length(x) > 0 && inherits(x[[1]], "expression"))
    #    stop("Internal error: lapply should call the function directly to prevent overhead.")
    return(lapply(x, FUN))
  }

  if (is.null(cl) && cores > 1) {
    if (.Platform$OS.type == "windows") {
      warning("'mclapply' is not available on Windows; reverting to single-core 'lapply'.")
      return(lapply(x, FUN))
    }
    return(parallel::mclapply(x, FUN, mc.cores = cores, mc.preschedule = preschedule))
  }

  if (inherits(cl, "cluster")) {
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
#' taken as input, and repeated indices of the first occurrence of each row
#' are returned.
#'
#' @returns A vector of row indices corresponding to the first ocurrence of a given row.
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
  if (length(di) == 0) return(1:nrow(m))
  mdup  <- m[di, , drop = FALSE]
  tmuniq <- t(m[-di, , drop = FALSE])
  mudup  <- unique(mdup)  # Unique duplicates

  ret <- rep(NA, nrow(m))
  ret[!d] <- 1:sum(!d)
  for (i in 1:nrow(mudup)) {
    muniq.match <- which(apply(tm == mudup[i, ], 2, all))
    ret[muniq.match[-1]] <- ret[muniq.match[1]]
  }
  return(ret)
}

