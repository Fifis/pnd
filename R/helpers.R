# There internal functions are not exported

# Check the number of cores, check if we are on CRAN testing; there is a policy of
# 2 cores only...
.checkCores <- function(cores) {
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
