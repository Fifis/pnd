# Some of these internal functions are not exported

safeF <- function(FUN, x, ...) tryCatch(FUN(x, ...), error = function(e) return(structure(NA, error = e)))

checkBadSafeF <- function(x) identical(as.logical(x), NA) && identical(names(attributes(x)), "error")

# Concatenate together with a comma between the terms
pasteAnd <- function(x) paste0(x, collapse = ", ")

# Print in scientific (exponential) format like 1.23e-03 for 0.001234
printE <- function(x, d = 2) sprintf(paste0("%1.", d, "e"), x)

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

