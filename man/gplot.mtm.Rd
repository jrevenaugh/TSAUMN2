% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtm_plots.R
\name{gplot.mtm}
\alias{gplot.mtm}
\title{Plot results from \code{\link{spec.mtm}}}
\usage{
gplot.mtm(
  mspec,
  jack = FALSE,
  period = F,
  trans = c("identity", "log10", "sqrt")
)
}
\arguments{
\item{mspec}{\code{\link{spec.mtm}} output}

\item{jack}{include jackknife confidence intervals if TRUE.  These are computed by \code{\link{spec.mtm}}
only if "jackknife = TRUE".  The confidence level defaults to 95%.  This can be changed in \code{\link{spec.mtm}} 
with the jkCIProb argument.}

\item{period}{display spectrum against period, rather than frequency, if TRUE}

\item{trans}{apply a transform to the spectral power; defaults to no transform}
}
\value{
a ggplot2 object
}
\description{
Attractive \link{ggplot2} plot of \code{\link{spec.mtm}} output
}
\examples{
s <- ts( sin( 2 * pi * 1:512 / 32 ) + rnorm( 512, sd = 0.25 ), deltat = 1 )
mspec <- spec.mtm( s, jackknife = TRUE, plot = FALSE )
mtm.plot( mspec, jack = TRUE, period = TRUE, trans = "log10" )
}
