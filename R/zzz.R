#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Parallel numerical derivatives v. 0.0.6 (2025-01-23).")

  # The number of cores is auto-detected based on the OS
  os <- Sys.info()[["sysname"]]
  if (os == "Linux") {
    cores <- sum(!duplicated(grep("^core id", readLines("/proc/cpuinfo"), value = TRUE)))
  } else if (os == "Darwin") {
    cores <- system("system_profiler SPHardwareDataType | grep -i 'Total Number of Cores'", TRUE)
    cores <- as.integer(gsub("[^0-9]", "", cores))
  } else {  # Unfortunately one cannot go to system files
    cores <- parallel::detectCores(logical = FALSE)
  }
  # Worst case: less than a quarter is returned
  if (is.null(cores) || is.na(cores)) cores <- max(1, floor(parallel::detectCores()/2) - 1)

  if (cores > 4) cores <- cores - 1  # Leaving some resources for the system
  msg <- paste0("Using up to ", cores, " cores for parallelism through mclapply forking on Linux.\n",
                "To use a custom cluster, create it via\n",
                "library(parallel)\ncl <- makeCluster(", cores, ")\n",
                "and export objects / load packages via clusterExport() and clusterEvalQ().")
  if (.Platform$OS.type != "unix") {
    msg <- "Using only one core (parallelism on Windows in development)."
  } else if (os == "Darwin") {
    msg <- gsub("Linux", "Mac", msg)
  }
  packageStartupMessage(msg)

  options(pnd.cores = getOption("pnd.cores", cores))
  options(pnd.preschedule = getOption("pnd.preschedule", TRUE))
  options(pnd.warn.vectorised = getOption("pnd.warn.vectorised", FALSE))
}
