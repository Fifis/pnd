#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Parallel numerical derivatives v. 0.0.3 (2023-06-01).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}
