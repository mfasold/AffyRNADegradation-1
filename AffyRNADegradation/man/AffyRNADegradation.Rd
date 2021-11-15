\name{RNADegradation}

\Rdversion{1.1}

\alias{RNADegradation}

\title{Compute degradation parameters and correct probe intensities}

\description{ Given an AffyBatch object, the function \code{RNADegradation}
  computes the probe positional bias and several statistical parameters,
  including a robust and accurate measure of RNA degradation. Probe
  intensities are corrected for the bias for each sample within the
  \code{AffyBatch} object.}

\usage{
RNADegradation(affyData,
               location.type = "index",
               location.file.dir = NULL,
               plot.images = FALSE)
}
\arguments{
  \item{affyData}{ an \code{AffyBatch} object containing the probe level microarray
    data.}
  \item{location.type}{ index-based probe alignment (x=k) if "index", or
    location-based alignment (x=L) if "absolute".}
  \item{location.file.dir}{ directory containing the probe location file(s).}
  \item{plot.images}{ if \code{TRUE} plots a set of debug images.}
}

% \details{
% \code{RNADegradation} is the principal function that determines the probe
% location bias for each sample, and applies the correction.
% }

\value{
  An \code{AffyDegradationBatch} object.
}

\author{Mario Fasold}

\seealso{AffyDegradationBatch}

\examples{
if (require(AmpAffyExample)) {
  # Get example data
  data(AmpData)

  affy.deg <- RNADegradation(AmpData[,4])
}
}

\keyword{methods}
