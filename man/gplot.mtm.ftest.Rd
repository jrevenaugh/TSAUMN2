% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtm_plots.R
\name{gplot.mtm.ftest}
\alias{gplot.mtm.ftest}
\title{Plot F-test results from \code{\link{spec.mtm}}}
\usage{
gplot.mtm.ftest(
  mspec,
  uci = 0.95,
  period = F,
  trans = c("identity", "log10", "sqrt")
)
}
\arguments{
\item{mspec}{\code{\link{spec.mtm}} output ("Ftest = TRUE" )}

\item{uci}{upper confidence interval to highlight in plot}

\item{period}{display F-test against period, rather than frequency, if TRUE}

\item{trans}{apply a transform to the spectral power; defaults to no transform}
}
\value{
a ggplot2 object
}
\description{
Attractive \link{ggplot2} plot of \code{\link{spec.mtm}} F-test output
}
\examples{
s <- ts( sin( 2 * pi * 1:512 / 32 ) + rnorm( 512 ), deltat = 1 )
mspec <- spec.mtm( s, Ftest = TRUE, plot = FALSE )
mtm.ftest( mspec )

}
