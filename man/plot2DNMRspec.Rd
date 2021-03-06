% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot2DNMRspec.R
\name{plot2DNMRspec}
\alias{plot2DNMRspec}
\title{Draw a 2D NMR Spectrum}
\usage{
plot2DNMRspec(peaks, x.range = c(0, 12), MHz = 300, ppHz = 1,
  type = "COSY", M = NULL, levels = seq(0.5, 1, by = 0.1), ...)
}
\arguments{
\item{peaks}{A data frame with the following columns: delta, mult
(multiplicity), J, area, pw.  Multiplicity should be given by a number, so
use 2 for a doublet.  J is in Hz (use 0 for singlets).  pw is the peak width
at half-height in Hz.}

\item{x.range}{A numeric vector of length 2 giving the ppm range desired.
Must be increasing.}

\item{MHz}{Integer.  The operating frequency of the instrument, in MHz.}

\item{ppHz}{Points per Hz: The number of data points per Hz to use in
calculating the spectrum (passed as argument \code{dd} to \code{makeSpec}).
The default (1) works well for 1H NMR spectra.
Note that this function uses Hz internally so that the \code{x.range}, which
is in ppm, is multiplied by \code{Mhz} before being sent to
\code{\link{makeSpec}}, and once there, \code{makeSpec} will multiply it by
\code{ppHz}.  Thus the total data points used is \code{ppHz * Mhz *
abs(diff(x.range))}.  This approach ensures that peaks are not distorted
when changing \code{x.range} for the same \code{peak.list}.}

\item{type}{The type of 2D spectrum desired.  One of \code{c("COSY", "TOCSY")}.}

\item{M}{An adjacency matrix indicating which peaks are coupled.
The order of rows and columns must be the same as in \code{peaks}.}

\item{levels}{A vector of levels for the contour plot.  Must be in (0...1).}

\item{\ldots}{Parameters to be passed to the plotting function.}
}
\value{
Returns a matrix.
}
\description{
This function simulates 2D NMR spectra.  Only 1st order coupling can be
handled -- there is currently no capacity for doublet of doublets and
other such peaks.  The field strength of the "instrument" is taken into
account.
}
\examples{

### ethyl 2-ethyl-3-oxobutyrate
### Set up data

peaks1 <- data.frame(
#             A     B     C     D     E     F
	delta = c(4.20, 3.34, 2.23, 1.88, 1.28, 0.94),
	mult = c(4, 3, 1, 5, 3, 3),
	J = c(14, 14, 0, 14, 14, 14),
	area = c(2, 1, 3, 2, 3, 3),
	pw = c(2, 2, 2, 2, 2, 2))

#              A, B, C, D, E, F
AM <- matrix(c(0, 0, 0, 0, 1, 0,  # A
               0, 0, 0, 1, 0, 0,  # B
               0, 0, 0, 0, 0, 0,  # C
               0, 1, 0, 0, 0, 1,  # D
               1, 0, 0, 0, 0, 0,  # E
               0, 0, 0, 1, 0, 0), # F
			   ncol = 6)

### 1D 1H NMR plot for reference
# CRAN checks will skip some examples to save time


jnk <- plotNMRspec(peaks = peaks1, x.range = c(0, 5), MHz = 500,
main = "1H NMR of ethyl 2-ethyl-3-oxobutyrate")

### 2D COSY plot

res <- plot2DNMRspec(peaks = peaks1, x.range = c(0, 5), MHz = 500, ppHz = 1, M = AM,
main = "COSY of ethyl 2-ethyl-3-oxobutyrate")

### 2D TOCSY plot

\dontrun{

res <- plot2DNMRspec(peaks = peaks1, x.range = c(0, 5), MHz = 500, ppHz = 1,
levels = c(0.85, 0.9, 0.95), type = "TOCSY",
main = "TOCSY of ethyl 2-ethyl-3-oxobutyrate")
}
}
\seealso{
\code{\link{makeSpec}}
}
\author{
Bryan A. Hanson, DePauw University. \email{hanson@depauw.edu}
}
\keyword{utilities}
