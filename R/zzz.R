#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Parallel numerical derivatives v. 0.0.4 (2024-06-10).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}
