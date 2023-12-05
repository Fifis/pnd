#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Parallel numerical gradient v. 0.0.2 (2023-12-06).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}

