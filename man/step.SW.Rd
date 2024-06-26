% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step-select.R
\name{step.SW}
\alias{step.SW}
\title{Stepleman--Winarsky automatic step selection}
\usage{
step.SW(
  FUN,
  x,
  h0 = 1e-05 * (abs(x) + (x == 0)),
  shrink.factor = 2,
  range = h0/c(1e+12, 1e-08),
  seq.tol = 1e-04,
  max.rel.error = .Machine$double.eps/2,
  maxit = 40L,
  diagnostics = FALSE,
  ...
)
}
\arguments{
\item{FUN}{Function for which the optimal numerical derivative step size is needed.}

\item{x}{Numeric scalar: the point at which the derivative is computed and the optimal step size is estimated.}

\item{h0}{Numeric scalar: initial step size, defaulting to a relative step of
slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).}

\item{shrink.factor}{A scalar greater than 1 that is used to divide the step size
during the search. The authors recommend 4, but this may be result in earlier
termination at slightly sub-optimal steps. Change to 3 for a faster search.}

\item{range}{Numeric vector of length 2 defining the valid search range for the step size.}

\item{seq.tol}{Numeric scalar: maximum relative difference between old and new
step sizes for declaring convergence.}

\item{max.rel.error}{Positive numeric scalar > 0 indicating the maximum relative
error of function evaluation. For highly accurate functions with all accurate bits
is equal to half of machine epsilon. For noisy functions (derivatives, integrals,
output of optimisation routines etc.), it is higher.}

\item{maxit}{Maximum number of algorithm iterations to avoid infinite loops.
Consider trying some smaller or larger initial step size \code{h0}
if this limit is reached.}

\item{diagnostics}{Logical: if \code{TRUE}, returns the full iteration history
including all function evaluations.}

\item{...}{Passed to FUN.}
}
\value{
A list similar to the one returned by \code{optim()}: \code{par} -- the optimal
step size found, \code{value} -- the estimated numerical first derivative (central
differences), \code{counts} -- the number of iterations (each iteration includes
four function evaluations), \code{abs.error} -- an estimate of the total
approximation error (sum of truncation and rounding errors),
\code{exitcode} -- an integer code indicating the termination status:
\code{0} indicates optimal termination within tolerance,
\code{2} is returned if there is no change in step size within tolerance,
\code{3} indicates a solution at the boundary of the allowed value range,
\code{4} signals that the maximum number of iterations was reached.
\code{message} is a summary message of the exit status.
If \code{diagnostics} is \code{TRUE}, \code{iterations} is a list
including the full step size search path,
argument grids, function values on those grids, estimated derivative values,
estimated error values, and monotonicity check results.
}
\description{
Stepleman--Winarsky automatic step selection
}
\details{
This function computes the optimal step size for central differences using the
\insertCite{stepleman1979adaptive}{pnd} algorithm.
}
\examples{
f <- function(x) x^4  # The derivative at 1 is 4
step.SW(x = 1, f)
step.SW(x = 1, f, h0 = 1e-9, diagnostics = TRUE) # Starting too low
# Starting somewhat high leads to too many preliminary iterations
step.SW(x = 1, f, h0 = 10, diagnostics = TRUE)
step.SW(x = 1, f, h0 = 1000, diagnostics = TRUE) # Starting absurdly high

f <- sin  # The derivative at pi/4 is sqrt(2)/2
step.SW(x = pi/4, f)
step.SW(x = pi/4, f, h0 = 1e-9, diagnostics = TRUE) # Starting too low
step.SW(x = pi/4, f, h0 = 0.1, diagnostics = TRUE) # Starting slightly high
# The following two example fail because the truncation error estimate is invalid
step.SW(x = pi/4, f, h0 = 10, diagnostics = TRUE)   # Warning
step.SW(x = pi/4, f, h0 = 1000, diagnostics = TRUE) # Warning
}
\references{
\insertAllCited{}
}
