\name{AffyRNADegradation}

\Rdversion{1.1}

\alias{AffyRNADegradation-package}
\alias{AffyRNADegradation}

\docType{package}

\title{AffyRNADegradation: analyze and correct RNA degradation effects}

\description{The \code{AffyRNADegradation} package helps in the assessment of RNA
  degradation effects in Affymetrix 3' expression arrays. The parameter d
  gives a robust and accurate measure of RNA integrity. The correction
  removes the probe positional bias, and thus improves comparability of
  samples that are affected by RNA degradation.}

\details{  
  \tabular{ll}{
Package: \tab AffyRNADegradation\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2011-10-13\cr
License: \tab GPL >=2\cr
}

The \code{RNADegradation} function performs both analysis and correction of
RNA degradation effects, returning an object of type
\code{AffyDegradationBatch}. The class provides accessor functions to obtain the
degradation parameter d and a \code{AffyBatch} object containing corrected probe
intensities. A more detailed look on the RNA degradation effects can be
gained through the tongs plot, the Dx plot and further statistics.
}

\author{Mario Fasold}

%\references{}

\examples{
if (require(AmpAffyExample)) {
  # Load example data (AmpData affybatch)
  data(AmpData)

  ## Compute and correct degradation for a single chip (for speed)
  affy.deg <- RNADegradation(AmpData[,4])

  ## Show degradation parameter d  
  d(affy.deg)

  ## Replace data with corrected data for further analysis
  AmpData <- afbatch(affy.deg)
}
}

\keyword{ package }
