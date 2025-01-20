newCluster <- function(cl = NULL, cores = 2) {
  if (cores == 1) return("lapply")
  windows <- .Platform$OS.type == "windows"
  max.cores <- getOption("pnd.cores", cores)
  if (cores > max.cores) warning(paste0("You requested more cores than you have physical ",
                                        "cores. Consider setting 'cores = ", max.cores, "'.\n"))
  if (!is.null(cl)) {
    if(!inherits(cl, "cluster"))  # Check if the cluster is valid
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

runParallel <- function(FUN, x, cores = 1, preschedule = FALSE, cl = NULL) {
  if (identical(cl, "lapply")) {
    ret <- lapply(x, FUN)
  } else if (identical(cl, "mclapply")) {
    ret <- parallel::mclapply(X = x, FUN = FUN, mc.silent = cores, mc.preschedule = preschedule)
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

