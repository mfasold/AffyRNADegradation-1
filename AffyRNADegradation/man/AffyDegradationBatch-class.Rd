\name{AffyDegradationBatch-class}

\Rdversion{1.1}

\docType{class}

\alias{AffyDegradationBatch-class}
\alias{AffyDegradationBatch,ANY}
\alias{AffyDegradationBatch}
\alias{d}
\alias{afbatch}
\alias{plotDx}
\alias{d,AffyDegradationBatch-method}
\alias{afbatch,AffyDegradationBatch-method}
\alias{plotDx,AffyDegradationBatch-method}

\title{Class AffyDegradationBatch}

\description{This class represents Affymetrix GeneChip probe level data
  that has been analysed and corrected for the probe location bias.}

\section{Objects from the Class}{
  Objects can be created using the function \code{RNADegradation}.
}

\section{Slots}{
  \describe{
    \item{\code{location.type}:}{Object of class \code{character}
  describing the type of probe position used for the analysis (probe index or probe location).}
    \item{\code{afbatch}:}{Object of class \code{AffyBatch} containing
  corrected probe level data.}
    \item{\code{stats}:}{Object of class \code{matrix} containg various
  statistical parameters from the analysis.}
    \item{\code{means.pm}:}{Object of class \code{matrix}
      containing the average PM probe intensites for probe locations of
      expressed genes.}
    \item{\code{means.mm}:}{Object of class \code{matrix}
      containing the average MM probe intensites for probe locations of
      expressed genes.}
  }
}

\section{Methods}{
  \describe{
   \item{afbatch}{\code{signature(x = "AffyDegradationBatch")}: returns
     as \code{AffyBatch} object containing corrected probe level data.}
   \item{d}{\code{signature(x = "AffyDegradationBatch")}: returns
     a vector with a measure of RNA degradation for each sample}
   \item{plotDx}{\code{signature(x = "AffyDegradationBatch")}: plots the
     probe location bias for all samples of the experiment.}   
}}

\author{Mario Fasold }

\examples{
if (require(AmpAffyExample)) {
  # Get example data
  data(AmpData)

  affy.deg <- RNADegradation(AmpData[,4])

  ## Plot degradation 
  plotDx(affy.deg)

  ## Get degradation measure
  d(affy.deg)

  ## Get AffyBatch object with corrected probe intensities
  afbatch(affy.deg)
}
}

\seealso{AffyBatch}

\keyword{classes}