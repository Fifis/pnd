% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{checkCores}
\alias{checkCores}
\title{Number of core checks and changes}
\usage{
checkCores(cores = NULL)
}
\arguments{
\item{cores}{Integer specifying the number of CPU cores used for parallel computation.
Recommended to be set to the number of physical cores on the machine minus one.}
}
\value{
An integer with the number of cores.
}
\description{
Number of core checks and changes
}
\examples{
checkCores()
checkCores(2)
suppressWarnings(checkCores(1000))
}
