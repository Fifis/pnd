% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gradient.R
\name{generateGrid}
\alias{generateGrid}
\title{Create a grid of points for a gradient / Jacobian}
\usage{
generateGrid(x, h, stencils, elementwise, vectorised)
}
\arguments{
\item{x}{Numeric vector or scalar: the point(s) at which the derivative is estimated.
\code{FUN(x)} must be finite.}

\item{h}{Numeric or character specifying the step size(s) for the numerical
difference or a method of automatic step determination (\code{"CR"}, \code{"CRm"},
\code{"DV"}, or \code{"SW"} to be used in \code{\link[=gradstep]{gradstep()}}). The default value is
described in \code{?GenD}.}

\item{stencils}{A list of outputs from \code{\link[=fdCoef]{fdCoef()}} for each coordinate of \code{x}.}

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
}
\value{
A list with points for evaluation, summation weights for derivative computation, and
indices for combining values.
}
\description{
Create a grid of points for a gradient / Jacobian
}
\examples{
generateGrid(1:4, h = 1e-5, elementwise = TRUE, vectorised = TRUE,
             stencils = lapply(1:4, function(a) fdCoef(acc.order = a)))
}
