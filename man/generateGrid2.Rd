% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hessian.R
\name{generateGrid2}
\alias{generateGrid2}
\title{Generate grid points for Hessians}
\usage{
generateGrid2(x, side, acc.order, h)
}
\arguments{
\item{x}{Numeric vector or scalar: the point(s) at which the derivative is estimated.
\code{FUN(x)} must be finite.}

\item{side}{Integer scalar or vector indicating the type of finite difference:
\code{0} for central, \code{1} for forward, and \code{-1} for backward differences.
Central differences are recommended unless computational cost is prohibitive.}

\item{acc.order}{Integer or vector of integers specifying the desired accuracy order
for each element of \code{x}.
The final error will be of the order \eqn{O(h^{\mathrm{acc.order}})}{O(h^acc.order)}.}

\item{h}{Numeric or character specifying the step size(s) for the numerical
difference or a method of automatic step determination (\code{"CR"}, \code{"CRm"},
\code{"DV"}, or \code{"SW"} to be used in \code{\link[=gradstep]{gradstep()}}). The default value is
described in \code{?GenD}.}
}
\value{
A list with elements:
\itemize{
\item \code{xlist}: a list of unique coordinate shifts,
\item \code{w}: the finite-difference weights (one per point),
\item \code{i1}, \code{i2}: integer vectors giving partial-derivative indices.
}
The length of each vector matches \code{xlist}.
}
\description{
Creates a list of unique evaluation points for second derivatives: both
diagonal (\eqn{\partial^2 / \partial x_i^2}{d^2 / dx^2}) and cross
(\eqn{\partial^2 / \partial x_i \partial x_j}{d^2 / dx_i dx_j}).
}
\examples{
generateGrid2(1:4, side = rep(0, 4), acc.order = c(2, 6, 4, 2),
              h = c(1e-5, 1e-4, 2e-5, 1e-6))
}
\seealso{
\code{\link[=GenD]{GenD()}}, \code{\link[=Hessian]{Hessian()}}.
}
