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
#' cl <- newCluster(cores = 1)
newCluster <- function(cl = NULL, cores = 2) {
  if (cores == 1) return("lapply")
  windows <- .Platform$OS.type == "windows"
  max.cores <- getOption("pnd.cores", cores)
  if (cores > max.cores) warning(paste0("You requested more cores than you have physical ",
                                        "cores. Consider setting 'cores = ", max.cores, "'.\n"))
  if (!is.null(cl)) {
    if (!inherits(cl, "cluster"))  # Check if the cluster is valid
      stop(paste0("The object passed as a cluster is not a cluster. ",
                  "Use 'cl <- makeCluster(", cores, ")' to create a proper one;",
                  "use 'clusterExport(objs)' and clusterEvalQ(library(lib)) to add user-loaded data/features."))
  } else if (windows) {
    cl <- parallel::makePSOCKcluster(cores)
    parallel::clusterExport(cl, setdiff(ls(envir = parent.env()), "cl"))
    core.pkgs <- c("parallel", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
    loaded.pkgs <- setdiff((.packages()), c(core.pkgs, "pnd"))
    if (l <- length(loaded.pkgs) > 0) {
      sapply(seq_len(l),
             function(i) eval(parse(text = paste0("clusterEvalQ(cl, library(", loaded.pkgs[i], "))"))))
    }
  } else {
    cl <- "mclapply"
  }
  return(cl)
}

#' Run a function in parallel over a list (internal use only)
#'
#' @param FUN A function of only one argument. If there are more arguments, use
#'   the \code{FUN2 <- do.call(FUN, c(list(x), ...))} annd call it.
#' @param x A list to parallelise the evaluation of \code{FUN} over.
#' @param cores Integer specifying the number of CPU cores used for parallel computation.
#'   Recommended to be set to the number of physical cores on the machine minus one.
#' @param preschedule Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
#'   or enables load balancing with \code{parLapplyLB()}. Recommended for functions that
#'   take less than 0.1 s per evaluation.
#' @param cl An optional user-supplied \code{cluster} object (created by \code{makeCluster}
#'   or similar functions). If not \code{NULL}, the code uses \code{parLapply()}
#'   (if \code{preschedule} is \code{TRUE}) or \code{parLapplyLB()} on that cluster on Windows
#'   and \code{mclapply} (fork cluster) on everything else.
#'
#' @returns
#' @export
#'
#' @examples
runParallel <- function(FUN, x, cores = 1, preschedule = FALSE, cl = NULL) {
  if (identical(cl, "lapply")) {
    ret <- lapply(x, FUN)
  } else if (identical(cl, "mclapply")) {
    ret <- parallel::mclapply(X = x, FUN = FUN, mc.cores = cores, mc.preschedule = preschedule)
  } else if (inherits(cl, "cluster")) {
    if (preschedule) {
      ret <- parallel::parLapply(cl = cl, X = x, fun = FUN)
    } else {
      ret <- parallel::parLapplyLB(cl = cl, X = x, fun = FUN)
    }
  } else {
    stop("'cl' should be 'lapply', 'mclapply', or a cluster object.")
  }
  return(ret)
}
