% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step-select.R, R/methods.R
\name{gradstep}
\alias{gradstep}
\alias{print.stepsize}
\alias{print.gradstep}
\title{Automatic step selection for gradients}
\usage{
gradstep(
  FUN,
  x,
  h0 = NULL,
  method = c("plugin", "SW", "CR", "CRm", "DV", "M", "K"),
  control = NULL,
  cores = 1,
  preschedule = getOption("pnd.preschedule", TRUE),
  cl = NULL,
  ...
)

\method{print}{stepsize}(x, ...)

\method{print}{gradstep}(x, ...)
}
\arguments{
\item{FUN}{Function for which the optimal numerical derivative step size is needed.}

\item{x}{Numeric vector or scalar: the point at which the derivative is computed
and the optimal step size is estimated.}

\item{h0}{Numeric vector or scalar: initial step size, defaulting to a relative step of
slightly greater than .Machine$double.eps^(1/3) (or absolute step if \code{x == 0}).}

\item{method}{Character indicating the method: \code{"CR"} for \insertCite{curtis1974choice}{pnd},
\code{"CRm"} for modified Curtis--Reid, "DV" for \insertCite{dumontet1977determination}{pnd},
\code{"SW"} \insertCite{stepleman1979adaptive}{pnd}, \code{"M"} for
\insertCite{mathur2012analytical}{pnd}, \code{"K"} for Kostyrka (2026, exerimental),
and \code{"plugin"} for the single-step plug-in estimator.}

\item{control}{A named list of tuning parameters for the method. If \code{NULL},
default values are used. See the documentation for the respective methods. Note that
full iteration history including all function evaluations is returned, but
different methods have slightly different diagnostic outputs.}

\item{cores}{Integer specifying the number of CPU cores used for parallel computation.
Recommended to be set to the number of physical cores on the machine minus one.}

\item{preschedule}{Logical: if \code{TRUE}, disables pre-scheduling for \code{mclapply()}
or enables load balancing with \code{parLapplyLB()}. Recommended for functions that
take less than 0.1 s per evaluation.}

\item{cl}{An optional user-supplied \code{cluster} object  (created by \code{makeCluster}
or similar functions). If not \code{NULL}, the code uses \code{parLapply()} (if \code{preschedule}
is \code{TRUE}) or \code{parLapplyLB()} on that cluster on Windows, and \code{mclapply}
(fork cluster) on everything else.}

\item{...}{Passed to FUN.}
}
\value{
A list similar to the one returned by \code{optim()} and made of
concatenated individual elements coordinate-wise lists: \code{par} -- the optimal
step sizes found, \code{value} -- the estimated numerical gradient,
\code{counts} -- the number of iterations for each coordinate,
\code{abs.error} -- an estimate of the total approximation error
(sum of truncation and rounding errors),
\code{exitcode} -- an integer code indicating the termination status:
\code{0} indicates optimal termination within tolerance,
\code{1} means that the truncation error (CR method) or the third derivative
(DV method) is zero and large step size is preferred,
\code{2} is returned if there is no change in step size within tolerance,
\code{3} indicates a solution at the boundary of the allowed value range,
\code{4} signals that the maximum number of iterations was reached.
\code{message} -- summary messages of the exit status.
\code{iterations} is a list of lists
including the full step size search path, argument grids, function values on
those grids, estimated error ratios, and estimated derivative values for
each coordinate.
}
\description{
Automatic step selection for gradients
}
\details{
We recommend using the Stepleman--Winarsky algorithm because it does not suffer
from over-estimation of the truncation error in the Curtis--Reid approach
and from sensitivity to near-zero third derivatives in the Dumontet--Vignes
approach. It really tries multiple step sizes and handles missing
values due to bad evaluations for inadequate step sizes really in a robust manner.
}
\examples{
gradstep(x = 1, FUN = sin, method = "CR")
gradstep(x = 1, FUN = sin, method = "CRm")
gradstep(x = 1, FUN = sin, method = "DV")
gradstep(x = 1, FUN = sin, method = "SW")
gradstep(x = 1, FUN = sin, method = "M")
gradstep(x = 1, FUN = sin, method = "K")
# Works for gradients
gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "CR")
gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "CRm")
gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "DV")
gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "SW")
gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "M")
gradstep(x = 1:4, FUN = function(x) sum(sin(x)), method = "K")
print(step.CR(x = 1, sin))
print(step.DV(x = 1, sin))
print(step.plugin(x = 1, sin))
print(step.SW(x = 1, sin))
print(step.M(x = 1, sin))
print(step.K(x = 1, sin))
f <- function(x) x[1]^3 + sin(x[2])*exp(x[3])
print(gradstep(x = c(2, pi/4, 0.5), f))
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link[=step.CR]{step.CR()}} for Curtis--Reid (1974) and its modification,
\code{\link[=step.plugin]{step.plugin()}} for the one-step plug-in solution,
\code{\link[=step.DV]{step.DV()}} for Dumontet--Vignes (1977),
\code{\link[=step.SW]{step.SW()}} for Stepleman--Winarsky (1979),
\code{\link[=step.M]{step.M()}} for Mathur (2012), and
\code{\link[=step.K]{step.K()}} for Kostyrka (2026).
}
