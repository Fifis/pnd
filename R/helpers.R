# There internal functions are not exported

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
  if (nzchar(chk) && chk == "TRUE") cores <- 2L
  return(cores)
}

# Try evaluating a function and return NA with an attribute in case of an error
.safeF <- function(FUN, x, ...) tryCatch(FUN(x, ...), error = function(e) return(structure(NA, error = "error")))

# Concatenate together with a comma between the terms
.pasteAnd <- function(x) paste0(x, collapse = ", ")

# Print in scientific (exponential) format like 1.23e-03 for 0.001234
.printE <- function(x, d = 2) sprintf(paste0("%1.", d, "e"), x)
