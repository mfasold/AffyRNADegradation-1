\name{GetTongs}

\Rdversion{1.1}

\alias{GetTongs}
\alias{PlotTongs}
\alias{PlotDegradationHooks}

\title{Generate and visualize the tongs plot}

\description{The tongs plot and the related degradation hook show the
  relationship between expression level and probe location bias.  They
  are required for the correct estimation of RNA degradation effects. }

\usage{
GetTongs(affyData, chip.idx = 1)
PlotTongs(tongs)
PlotDegradationHooks(affyData, ...)
}
                            
\arguments{
  \item{affyData}{ an \code{AffyBatch} object.}
  \item{chip.idx}{ index of the sample to compute the tongs for.}
  \item{tongs}{ the tongs plot data table.}
  \item{...}{ optional graphical parameters.}
}

\value{
  \item{Tongs}{ a table containing Tongs plot values ordered by
    expression level.}
}

\author{Mario Fasold}

\examples{
if (require(AmpAffyExample)) {
  # Get example data
  data(AmpData)

  tongs <- GetTongs(AmpData, chip.idx = 4)
  PlotTongs(tongs)

  PlotDegradationHooks(AmpData[,c(3,4)])
}
}

\keyword{ methods }
