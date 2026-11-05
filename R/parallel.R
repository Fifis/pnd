# Only clusters created by pnd are retained here.
.pndCluster <- new.env(parent = emptyenv())
.pndCluster$cl <- NULL
.pndCluster$pid <- NULL

getParallelCluster <- function(cores, cl = NULL) {
  if (is.null(cl)) cl <- parallel::getDefaultCluster()
  if (!is.null(cl)) {
    if (!inherits(cl, "cluster"))
      stop("The object passed as a cluster is not a cluster. Use 'parallel::makePSOCKcluster()'.")
    return(cl)
  }
  if (cores == 1 || isTRUE(getOption("pnd.parallel.worker", FALSE))) return(NULL)

  backend <- match.arg(getOption("pnd.parallel", "auto"), c("auto", "PSOCK", "mclapply"))
  if (backend == "auto")
    backend <- if (.Platform$OS.type == "windows") "PSOCK" else "mclapply"
  if (backend == "mclapply") {
    if (.Platform$OS.type == "windows") stop("'mclapply' is unavailable on Windows; use 'PSOCK'.")
    return(NULL)
  }

  # Forked children must not use or stop their parent's socket connections.
  if (!identical(.pndCluster$pid, Sys.getpid())) .pndCluster$cl <- NULL
  if (is.null(.pndCluster$cl) || length(.pndCluster$cl) < cores) {
    stopParallel()
    cl <- parallel::makePSOCKcluster(cores)
    # Workers need the same library paths, including temporary check libraries.
    setup <- function(paths) {
      .libPaths(paths)
      # Nested derivatives use these workers locally unless given a separate cluster.
      options(pnd.parallel.worker = TRUE)
      NULL
    }
    environment(setup) <- baseenv()
    tryCatch(parallel::clusterCall(cl, setup, .libPaths()), error = function(e) {
      try(parallel::stopCluster(cl), silent = TRUE)
      stop(e)
    })
    .pndCluster$cl <- cl
    .pndCluster$pid <- Sys.getpid()
  }
  # Smaller requests reuse a subset, retaining the full pool for later calls.
  .pndCluster$cl[seq_len(cores)]
}

#' Stop automatically created parallel workers
#'
#' Stops the reusable PSOCK cluster created by pnd. Supplied clusters and clusters
#' registered through \code{parallel::setDefaultCluster()} are left to the caller.
#' The next parallel call creates workers again if needed. See [runParallel()].
#'
#' @returns \code{NULL}, invisibly.
#' @export
#' @examples
#' stopParallel()
stopParallel <- function() {
  cl <- .pndCluster$cl
  .pndCluster$cl <- NULL
  if (!is.null(cl) && identical(.pndCluster$pid, Sys.getpid())) parallel::stopCluster(cl)
  invisible(NULL)
}
