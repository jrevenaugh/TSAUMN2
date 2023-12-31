% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spec.ls.R
\name{spec.ls}
\alias{spec.ls}
\title{Compute Lomb-Scargle periodogram of (un)evenly spaced data}
\usage{
spec.ls(
  x,
  t = NULL,
  over = 1,
  demean = TRUE,
  detrend = FALSE,
  plot = TRUE,
  ...
)
}
\arguments{
\item{x}{The data to be analysed. x can be either a two-column numerical dataframe or matrix, with sampling times in columnn 1 and measurements in column 2, a single numerical vector containing measurements, or a single vector ts object (which will be converted to a numerical vector).}

\item{t}{If x is a single vector, t can be provided as a numerical vector of equal length containing sampling times. If x is a vector and times is NULL, the data are assumed to be equally sampled and times is set to 1:length(x).}

\item{over}{The oversampling factor. Must be an integer >= 1. Larger values of over lead to finer scanning of frequencies.}

\item{demean}{remove mean from timeseries prior to spectral estimation}

\item{detrend}{remove linear trend from timeseries prior to spectral estimation}

\item{plot}{Logical.  If plot = TRUE, the spectrum is plotted.}

\item{...}{Additional arguments passed to gplot.mtm.}
}
\value{
object of class spec with the following list items:
  \item{"freq"}{A vector with spectrum frequencies}
  \item{"spec"}{A vector with spectral power estimates corresponding to "freq"}
  \item{"series"}{Name of input time series}
  \item{"method"}{Method name: "Lomb-Scargle"}
}
\description{
Computes the Lomb-Scargle periodogram for a time series with irregular (or regular) sampling intervals.
}
\examples{
x <- rnorm( 256 )
s <- spec.ls( x, 1:256, plot = TRUE )
}
