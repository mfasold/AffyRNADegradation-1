## Creates tongs table
ComputeTongs <- function(x, intensity, sigma, probeset.id, min.x = 1, max.x = 11, mvg.avg.size = 450) {
  ## Compute boundaries (x) of upper and lower tongs (=upper and lower third of x)
  tongs.probes.per.probeset <- floor((max.x - min.x) / 3)
  p3.probes.max.x <- min.x + tongs.probes.per.probeset - 1 # 3' probes
  p5.probes.max.x <- max.x - tongs.probes.per.probeset + 1 # 5' probes

  ## Remove all probes outside min.x and max.x boundaries
  subset.x <- x >= min.x & x <= max.x 
  x <- x[subset.x]
  intensity <-  intensity[subset.x]
  sigma <- sigma[subset.x]
  probeset.id <-  probeset.id[subset.x]
   
  ## Firstly, select probesets where upper and lower tongs can be computed, that is, there exists 
  ## a probe within the upper and lower third of x available
  min.x.ps <- tapply(x, probeset.id, min)
  max.x.ps <- tapply(x, probeset.id, max)
  tongs.is.possible <- min.x.ps <= p3.probes.max.x & max.x.ps >= p5.probes.max.x
  tongs.is.possible.ps <- rownames(tongs.is.possible)[tongs.is.possible == T]

  subset.x <- probeset.id %in% tongs.is.possible.ps
  x <- x[subset.x]
  intensity <-  intensity[subset.x]
  sigma <- sigma[subset.x]
  probeset.id <-  probeset.id[subset.x]

  intensity <- log10(intensity) # Use log 10

  ## Compute probset-based averages
  sigma.i <- tapply(sigma, probeset.id, ProbesetAverage) 
  # mean von gleichen werten ist der wert... (falls sigma schon probeset average ist)
  average.i <- tapply(intensity, probeset.id, ProbesetAverage)
  is.p3.x <- x <= p3.probes.max.x 
  p3.ps.i <- tapply(intensity[is.p3.x], probeset.id[is.p3.x], ProbesetAverage)
  is.p5.x <- x >= p5.probes.max.x
  p5.ps.i <- tapply(intensity[is.p5.x], probeset.id[is.p5.x], ProbesetAverage)

  mean.x <- tapply(x, probeset.id, mean)
  
  tongs <- data.frame(sigma = sigma.i, mean.x = mean.x,
    p3.i = p3.ps.i, p5.i = p5.ps.i, average.i = average.i, row.names = rownames(mean.x))

  # Apply moving average values
  tongs <- tongs[order(tongs$sigma),]
  tongs[,"p5.tongs"] <- filter(tongs[,"p5.i"] - tongs[,"average.i"], rep(1/mvg.avg.size,mvg.avg.size))
  tongs[,"p3.tongs"] <- filter(tongs[,"p3.i"] - tongs[,"average.i"], rep(1/mvg.avg.size,mvg.avg.size))
  tongs[,"hook2"] <- filter(tongs[,"p3.i"] - tongs[,"p5.i"], rep(1/mvg.avg.size,mvg.avg.size))
  
  tongs <- na.omit(tongs) # remove NA rows created by filter
  tongs[,"hook"] <- tongs[,"p3.tongs"] - tongs[,"p5.tongs"]

  return(tongs)
}

EstimateHookParams <- function(tongs) {
  tongs.loess <- loess(hook ~ sigma, tongs, span=2/3)
  
  hook.params <- list()
  hook.params$max.raw <- max(tongs$hook) 
  hook.params$max.fit <- max(tongs.loess$fit)                                        

  tongs$hook <- tongs.loess$fit / max(tongs.loess$fit)
  tongs$hook <- pmax(tongs$hook, 0) # forbid negatives

  hook.params$hook.max.sigma <- tongs[min(which(tongs$hook == max(tongs$hook))),"sigma"]

  ## Estimate non-specific threshold as maximum sigma where (smoothed) hook is smaller than 0.02
  hook.params$ns.threshold <- max(tongs$sigma[ tongs.loess$fit < 0.02 & tongs$sigma < hook.params$hook.max.sigma])

  ## Alternate estimate: simply consider 30% smallest EMs as non-specific
  hook.params$ns.threshold2 <- quantile(tongs$sigma, 0.30) 

  # Make sure ns.threshold parameter is available
  if (is.nan(hook.params$ns.threshold) || is.infinite(hook.params$ns.threshold))
    hook.params$ns.threshold <- hook.params$ns.threshold2

  return(hook.params)
}

# Currently only give tongs for index based variant, x=k
GetTongs <- function(affyData, chip.idx = 1) {
  probe.info <- getProbeInfo(affyData, "index")

  intensity <- as.vector(pm(affyData[,chip.idx]))
  sigma     <- ave(intensity, list(probe.info$probeset.id), FUN=ProbesetLogAverage)
  
  ## Create tongs and estimate parameters
  tongs <- ComputeTongs(probe.info$x, intensity, sigma, probe.info$probeset.id, 1, 11)
    
  ##return(new("Tongs", tongs));
  return(tongs)
}

# Plotting helper
PlotTongs <- function(tongs) {
 ymax <- 0.35
 ymin <- -0.35 
 plot(tongs[,"sigma"], tongs[,"p5.tongs"], type="l", ylim=c(ymin, ymax),
      xlab=expression(Sigma), ylab=expression(Sigma[subset] - Sigma),
      main=paste("Tongs plot"), panel.first = grid())
 lines(tongs[,"sigma"], tongs[,"p3.tongs"], col=2)
 legend("topleft", c("5' subset", "3' subset"), fill=seq(1,2))
}

# Plot degradatio hook
# (difficult to provide sensible defaults for axis limits)
PlotDegradationHooks <- function(affyData, ...) {
  colors = rainbow(length(affyData))
  plot(0, xlim=c(1.4, 4.1), ylim=c(-0.03, 0.5), type="n", 
       xlab=expression(Sigma), ylab=expression(Delta), ...)

  # Compute and plot degradation hook separately for each sample
  for (chip.idx in seq_along(affyData)) {
    tongs <- GetTongs(affyData, chip.idx=chip.idx)
    lines(tongs$sigma, tongs$hook, col=colors[chip.idx], lwd=2, ...)
  }
  legend("topleft", sampleNames(affyData), fill=colors, ...)
}
