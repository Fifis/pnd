#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Parallel numerical derivatives v. 0.0.5 (2025-01-15).")
  packageStartupMessage("This is a pre-release version. Core functions subject to change.")
}
