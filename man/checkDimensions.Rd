% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gradient.R, R/methods.R
\name{checkDimensions}
\alias{checkDimensions}
\alias{print.checkDimensions}
\title{Determine function dimensionality and vectorisation}
\usage{
checkDimensions(
  FUN,
  x,
  f0 = NULL,
  func = NULL,
  elementwise = NA,
  vectorised = NA,
  multivalued = NA,
  deriv.order = 1,
  acc.order = 2,
  side = 0,
  h = NULL,
  zero.tol = sqrt(.Machine$double.eps),
  cores = 1,
  preschedule = TRUE,
  cl = NULL,
  ...
)

\method{print}{checkDimensions}(x, ...)
}
\arguments{
\item{FUN}{A function returning a numeric scalar or a vector whose derivatives are to be
computed. If the function returns a vector, the output will be a Jacobian.}

\item{x}{Numeric vector or scalar: the point(s) at which the derivative is estimated.
\code{FUN(x)} must be finite.}

\item{f0}{Optional numeric: if provided, used to determine the vectorisation type
to save time. If FUN(x) must be evaluated (e.g. second derivatives), saves one evaluation.}

\item{func}{For compatibility with \code{numDeriv::grad()} only. If instead of
\code{FUN}, \code{func} is used, it will be reassigned to \code{FUN} with a warning.}

\item{elementwise}{Logical: is the domain effectively 1D, i.e. is this a mapping
\eqn{\mathbb{R} \mapsto \mathbb{R}}{R -> R} or
\eqn{\mathbb{R}^n \mapsto \mathbb{R}^n}{R^n -> R^n}. If \code{NA},
compares the output length ot the input length.}

\item{vectorised}{Logical: if \code{TRUE}, the function
is assumed to be vectorised: it will accept a vector of parameters and return
a vector of values of the same length. Use \code{FALSE} or \code{"no"}  for
functions that take vector arguments and return outputs of arbitrary length (for
\eqn{\mathbb{R}^n \mapsto \mathbb{R}^k}{R^n -> R^k} functions). If \code{NA},
checks the output length and assumes vectorisation if it matches the input length;
this check is necessary and potentially slow.}

\item{multivalued}{Logical: if \code{TRUE}, the function is assumed to return vectors longer
than 1. Use \code{FALSE} for element-wise functions. If \code{NA}, attempts inferring it from
the function output.}

\item{deriv.order}{Integer or vector of integers indicating the desired derivative order,
\eqn{\mathrm{d}^m / \mathrm{d}x^m}{d^m/dx^m}, for each element of \code{x}.}

\item{acc.order}{Integer or vector of integers specifying the desired accuracy order
for each element of \code{x}.
The final error will be of the order \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.}

\item{side}{Integer scalar or vector indicating the type of finite difference:
\code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
Central differences are recommended unless computational cost is prohibitive.}

\item{h}{Numeric or character specifying the step size(s) for the numerical
difference or a method of automatic step determination (\code{"CR"}, \code{"CRm"},
\code{"DV"}, or \code{"SW"} to be used in \code{\link[=gradstep]{gradstep()}}). The default value is
described in \code{?GenD}.}

\item{zero.tol}{Small positive integer: if \code{abs(x) >= zero.tol}, then, the automatically
guessed step size is relative (\code{x} multiplied by the step), unless an auto-selection
procedure is requested; otherwise, it is absolute.}

\item{cores}{Integer specifying the number of CPU cores used for parallel computation.
Recommended to be set to the number of physical cores on the machine minus one.}

\item{preschedule}{Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
or enables load balancing with \code{parLapplyLB()}. Recommended for functions that
take less than 0.1 s per evaluation.}

\item{cl}{An optional user-supplied \code{cluster} object  (created by \code{makeCluster}
or similar functions). If not \code{NULL}, the code uses \code{parLapply()} (if \code{preschedule}
is \code{TRUE}) or \code{parLapplyLB()} on that cluster on Windows, and \code{mclapply}
(fork cluster) on everything else.}

\item{...}{Ignored.}
}
\value{
A named logical vector indicating if a function is element-wise or not,
vectorised or not, and multivalued or not.
}
\description{
Determine function dimensionality and vectorisation
}
\details{
The following combinations of parameters are allowed, suggesting specific input and
output handling by other functions:\tabular{lll}{
    \tab \code{elementwise} \tab \code{!elementwise} \cr
   \code{!multivalued}, \code{vectorised} \tab \code{FUN(xgrid)} \tab \emph{(undefined)} \cr
   \code{!multivalued}, \code{!vectorised} \tab \verb{[mc]lapply(xgrid, FUN)} \tab \verb{[mc]lapply} gradient \cr
   \code{multivalued}, \code{vectorised} \tab \emph{(undefined)} \tab \code{FUN(xgrid)} Jacobian \cr
   \code{multivalued}, \code{!vectorised} \tab \emph{(undefined)} \tab \verb{[mc]lapply} Jacobian \cr
}


Some combinations are impossible: multi-valued functions cannot be element-wise,
and single-valued vectorised functions must element-wise.

In brief, testing the input and output length and vectorisation capabilities may result in five
cases, unlike 3 in \code{numDeriv::grad()} that does not provide checks for Jacobians.
}
\examples{
checkDimensions(sin, x = 1:4, h = 1e-5)  # Rn -> Rn vectorised
checkDimensions(function(x) integrate(sin, 0, x)$value, x = 1:4, h = 1e-5)  # non vec
checkDimensions(function(x) sum(sin(x)), x = 1:4, h = 1e-5)  # Rn -> R, gradient
checkDimensions(function(x) c(sin(x), cos(x)), x = 1, h = 1e-5)  # R -> Rn, Jacobian
checkDimensions(function(x) c(sin(x), cos(x)), x = 1:4, h = 1e-5)  # vec Jac
checkDimensions(function(x) c(integrate(sin, 0, x)$value, integrate(sin, -x, 0)$value),
                 x = 1:4, h = 1e-5)  # non-vectorised Jacobian
}
