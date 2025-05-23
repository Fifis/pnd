#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Parallel numerical derivatives v.0.1.0 (2025-05-20).")

  # The number of cores is auto-detected based on the OS
  os <- Sys.info()[["sysname"]]
  if (os == "Linux") {
    cores <- sum(!duplicated(grep("^core id", readLines("/proc/cpuinfo"), value = TRUE)))
  } else if (os == "Darwin") {
    # Accepting only unambiguous values
    cores <- tryCatch(as.integer(system("/usr/sbin/sysctl -n hw.physicalcpu", TRUE)),
                      warning = function(e) return(NULL),
                      error = function(e) return(NULL))
  } else {  # Unfortunately one cannot dig inside system files
    cores <- parallel::detectCores(logical = FALSE)
  }
  # Worst case: less than a quarter is returned
  if (is.null(cores) || is.na(cores)) cores <- max(1, floor(parallel::detectCores()/2) - 1)

  msg <- paste0(cores, " physical cores for parallelism through mclapply forking are available on Linux.\n")
  if (os == "Windows") msg <- paste0(msg, "Create and register a default cluster first.")
  if (.Platform$OS.type != "unix") {
    msg <- gsub("mclapply forking are available on Linux", "PSOCK cluster workers are available on Windows", msg)
  } else if (os == "Darwin") {
    msg <- gsub("Linux", "Mac", msg)
  }
  packageStartupMessage(msg)

  max.cores <- parallel::detectCores()  # Max. useful # of cores assuming availability

  # Providing global options for faster invocation later
  options(pnd.cores = getOption("pnd.cores", cores))
  options(pnd.max.cores = getOption("pnd.max.cores", max.cores))
  options(pnd.preschedule = getOption("pnd.preschedule", TRUE))
  options(pnd.warn.vectorised = getOption("pnd.warn.vectorised", FALSE))
}
