# Some of these internal functions are not exported

safeF <- function(FUN, x, ...) tryCatch(FUN(x, ...), error = function(e) return(structure(NA, error = "error")))

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
#' checkCores(100)
checkCores <- function(cores = NULL) {
  if (is.null(cores)) cores <- max(floor(parallel::detectCores()/2 - 1), 1L)

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")  # Limit to 2 cores for CRAN checks
  if (nzchar(chk) && chk == "TRUE") cores <- min(cores, 2L)
  return(cores)
}

# Constructs an expression local({z <- x[i]; FUN(z, a = a, b = b, ...) }, envir = e)
# Borrowed from optimParallel
getExpr <- function(i, fun, dots, fname = "FUN") {
  ex <- paste0(fname, "(z")
  if (!is.primitive(fun)) {
    ff <- formals(fun)
    if (names(ff)[1] != "...") ff <- ff[-1]
    if (all(names(ff) != "..."))
      ff <- ff[names(ff) %in% names(dots)]
    else
      ff <- dots
    if(length(ff) >= 1){
      ex <- paste0(ex, ",")
      moreArgs <- paste(lapply(names(ff), function(x) paste0(x, "=", x)), collapse = ", ")
      ex <- paste0(ex, moreArgs)
    }
  }
  fc <- paste0(ex, ")")

  # Refer to .GlobalEnv$e so that eval() find the exported environment
  ex <- paste0("local({z <- x[", i, "]; ", fc, "}, envir = .GlobalEnv$e)")
  parse(text = ex)
}


#' Create a cluster inside pnd-specific functions (internal use only)
#'
#' @param cl \code{NULL} = return a PSOCK cluster with exported objects on Windows
#'   or the \code{'lapply' / 'mclapply'} string on everything else
#' @param cores If 1, return \code{'lapply'}; if more, return a cluster on Windows
#'   or \code{'mclapply'} on everything else.
#'
#' @returns \code{'lapply'}, \code{'mclapply'}, or a cluster object.
#' @export
#'
#' @examples
#' cl <- checkOrCreateCluster(cores = 1)  # "lapply"
#' cl <- tryCatch(checkOrCreateCluster("mclapply 2"), error = function(e) return(NULL))
#' print(cl)
checkOrCreateCluster <- function(cl = NULL, cores = 1) {
  if (is.character(cl) && (cl != "lapply") && !grepl("^mclapply \\d+$", cl))
    stop(paste0("The string for the cluster must be either 'lapply' or 'mclapply N' (N being ",
                "the desired number of cores). On Windows, must be a valid cluster."))
  if (is.character(cl) && (cl == "mclapply 1")) return("lapply")
  if (is.character(cl) && cl == "lapply") return("lapply")
  # If it is a valid cluster and not a string, return it
  if (!is.null(cl) && inherits(cl, "cluster")) return(cl)
  # If 1 core was requested, use lapply
  if ((!is.null(cores)) && cores == 1) return("lapply")

  windows <- .Platform$OS.type == "windows"

  if (is.character(cl) && grepl("^mclapply ", cl)) {
    s <- strsplit(cl, " ")[[1]]
    cores <- as.integer(s[2]) # Ignore cores because the cluster is more important
    if (windows) {  # Skip returning the fork request
      if (cores > 1) warning("'mclapply' unavailable on Windows; use 'makeCluster()' instead.")
      return("lapply")
    } else {
      return(cl)
    }
  }

  # At this point, not a string and not a cluster implies rubbish input
  if (!is.null(cl))
    stop(paste0("The object passed as a cluster is not a cluster. ",
                "Use 'cl <- makeCluster(", cores, ")' to create a proper one."))

  cl <- parallel::getDefaultCluster() # Maybe a default cluster exists -- return it
  if (!is.null(cl)) return(cl)

  max.cores <- getOption("pnd.cores", cores)
  if (cores > max.cores) warning(paste0("You requested more cores than you have physical ",
                                        "cores. Consider setting 'cores = ", max.cores, "'.\n"))

  # Finally, if cl is still NULL, create it or schedule forking
  if (windows) {
    cl <- parallel::makePSOCKcluster(cores)
    parallel::setDefaultCluster(cl)
    message("No cluster specified; created a PSOCK cluster with ", cores, " workers on Windows.")
  } else {
    cl <- paste("mclapply", cores)
  }

  return(cl)
}

setupParallelEnv <- function(FUN, cl, ...) {
  dots <- list(...)
  e <- list2env(dots)
  FUNe <- FUN
  environment(FUNe) <- e
  assign("FUN", FUNe, envir = e)
  if (inherits(cl, "cluster")) {
    parallel::clusterExport(cl, "e", envir = list2env(list(e = e)))
  } else {
    assign("e", e, envir = .GlobalEnv)
  }
  list(FUN = FUNe, e = e, dots = dots)
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
#' @param cl A string \code{"lapply"}, \code{"mclapply X"} (where \code{X} is the number
#'   of cores, e.g. \code{"mclapply 8"}), or an optional user-supplied \code{cluster} object
#'   (created by \code{makeCluster} or similar functions). If not \code{NULL},
#'   the code uses \code{parLapply()} (if \code{preschedule} is \code{TRUE}) or
#'   \code{parLapplyLB()} on that cluster on Windows, and \code{mclapply}
#'   (fork cluster) on everything else.
#'
#'
#' @returns The value that `lapply(x, FUN)` would have returned.
#' @export
#'
#' @examples
#' fslow <- function(x) Sys.sleep(x)
#' x <- rep(0.05, 6)
#' cl <- checkOrCreateCluster(cores = 2)
#' print(t1 <- system.time(runParallel(fslow, x, cl = "lapply")))
#' print(t2 <- system.time(runParallel(fslow, x, cl = cl)))
#' cat("Parallel overhead at 2 cores: ", round(t2[3]*200/t1[3]-100), "%\n", sep = "")
#' if (inherits(cl, "cluster")) parallel::stopCluster(cl)
runParallel <- function(FUN, x, preschedule = FALSE, cl = NULL) {
  if (identical(cl, "lapply") || is.null(cl)) {
    if (is.list(x) && length(x) > 0 && inherits(x[[1]], "expression"))
      stop("Internal error: lapply should call the function directly to prevent overhead.")
    return(lapply(x, FUN))
  } else if (is.character(cl) && grepl("^mclapply\\s+\\d+$", cl)) {
    if (is.list(x) && length(x) > 0 && inherits(x[[1]], "expression"))
      stop("Internal error: mclapply should call the function directly to prevent overhead.")
    cores <- as.integer(strsplit(cl, " ")[[1]][2])
    return(parallel::mclapply(x, FUN, mc.cores = cores, mc.preschedule = preschedule))
  } else if (inherits(cl, "cluster")) {
    if (is.list(x) && length(x) > 0 && inherits(x[[1]], "expression")) {
      # If the first element of x is an expression, assume that every element is.
      if (preschedule)
        return(parallel::parLapply(cl, x, eval))
      else
        return(parallel::parLapplyLB(cl, x, eval))
    } else {  # x contains values for direct application
      if (preschedule)
        return(parallel::parLapply(cl, x, FUN))
      else
        return(parallel::parLapplyLB(cl, x, FUN))
    }
  }
  stop("'cl' should be 'lapply', 'mclapply', or a cluster object.")
}

