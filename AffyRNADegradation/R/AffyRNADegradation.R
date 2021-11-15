RNADegradation <- function(affyData,
                           location.type = "index",
                           location.file.dir = NULL,
                           plot.images = FALSE) { 
  ## Get probe Info
  probe.info <- getProbeInfo(affyData, location.type, location.file.dir)

  if (location.type == "index") { ##
    scale.x  <- 1 
    min.x    <- 1
    max.x    <- 11
    n.breaks <- 11
    positions <- min.x:n.breaks
  } else {
    scale.x  <- 60 # enables to use the same defaults for fitting in absolute and index
    min.x    <- 0
    max.x    <- 600
    n.breaks <- 12
    positions <- seq((max.x-min.x)/n.breaks/2, by = (max.x-min.x)/n.breaks, length.out = n.breaks)

    ## Print some statistics about absolute probe locations
    message(round(sum(!is.na(probe.info$x )) / nrow(probe.info), 2)*100,
        "% of all probes can be aligned to their target mRNA",
        " for the current chip type (", cdfName(affyData), ").", sep="")
  }

  ## Get basic probe and probeset information
  x <- probe.info$x
  probeset.id   <- probe.info$probeset.id
  probeset.size <- ave(probe.info$probeset.id, list(probe.info$probeset.id), FUN=length)

  ## Set up statistcs table  
  stats.all <- data.frame(d.inf = rep(0, length(affyData)),  lambda.inv  = 0,  shift.y  = 0,
                          tongs.max.raw = 0, tongs.max.fit = 0, hook.max.sigma = 0, 
                          ns.threshold = 0, placeholder = 0, decay = 0,
                          row.names = sampleNames(affyData))
  means.pm.all <- matrix(0, n.breaks, length(affyData), dimnames = list(x.center = positions, sample = sampleNames(affyData)))
  means.mm.all <- matrix(0, n.breaks, length(affyData), dimnames = list(x.center = positions, sample = sampleNames(affyData)))


  ## Iterate over chips
  for(chip.idx in 1:length(affyData)) {
    currentSample = sampleNames(affyData)[chip.idx]
    message("Computing degradation for chip", chip.idx, ":", currentSample,"\n")

    ## Create main table of chip probes
    intensity <- as.vector(pm(affyData[,chip.idx]))
    sigma     <- ave(intensity, list(probe.info$probeset.id), FUN=ProbesetLogAverage)

    ## Create tongs and estimate parameters
    tongs <- ComputeTongs(x, intensity, sigma, probeset.id, min.x, max.x)
    hook.params <- EstimateHookParams(tongs)
    if (plot.images) PlotTongs(tongs)

    ## Compute d(x) once (Means and Fit)
    h.max <- hook.params$hook.max.sigma
    subset.probes <- (probeset.size == 11) & (h.max - 0.4 < sigma) & (sigma < h.max + 0.2)
    means <- getGeneralDegradationMeans(x[subset.probes], intensity[subset.probes], min.x, max.x)
    d.x <- means$i / means$i[1]
    theoretic.d.s <- TryToFitDecayFunction((means$x - min.x)/scale.x, d.x, d.x[length(d.x)], 0.6) 

    decay = (d.x[length(d.x)] + d.x[length(d.x)-1]) / 2
    
    ## Compute scaling function for chip
    x.limited = pmin(x, max.x) # probes located more far away than max.x are corrected using max.x
    theoretic.d.s.all  <- with(theoretic.d.s,  DecayFunction((x.limited - min.x)/scale.x, d1, d2, shiftY))
    tongs.loess <- loess(hook ~ sigma, tongs, span=2/3)
    
    ## Define scaling function such that the average scaling of the region that has been used to
    ## estimate D(s) is 1. 
    scale.f.y <- pmax(0, tongs.loess$fit) # scaling must never be smaller than 0
    average.scaling.in.specific <- mean(scale.f.y[h.max - 0.4 < tongs$sigma & tongs$sigma < h.max + 0.2])
    scale.f.y <- scale.f.y / average.scaling.in.specific

    ## Special handling for PrimeView arrays which do not contain mismatch probes
    if (cdfName(affyData) == "PrimeView") {
      probe.types <- c("pm")      
    } else {
      probe.types <- c("mm","pm")
    }

    ## Separately correct pm and mm      
    for(probe.type in probe.types) { # mm before pm, to avoid using previously corrected pm values

      ## Assign PM or MM accessor functions, for convenience
      if (probe.type == "pm") {
        pmOrMm <- pm
        `pmOrMm<-` <- `pm<-`
      } else {
        pmOrMm <- mm
        `pmOrMm<-` <- `mm<-`
      }
      intensity <- as.vector(pmOrMm(affyData[,chip.idx]))

      ## Use separate sigma for PM AND MM
      sigma <- ave(intensity, list(probeset.id), FUN=ProbesetLogAverage)

      means <- getGeneralDegradationMeans(x[subset.probes], intensity[subset.probes], min.x, max.x)
      # print(means)
      if (plot.images) PlotDx(means, min.x, max.x, scale.x)

      ## Save PM/MM mean statistics
      if (probe.type == "pm") {
        means.pm.all[,chip.idx] <- means$i 
      } else {
        means.mm.all[,chip.idx] <- means$i
      }

      ## Apply scaling function to sigma and x values
      if (plot.images) plot(tongs$sigma, scale.f.y, main=paste("Scaling function, ", probe.type))
      scaling.x <- approx(tongs$sigma, scale.f.y, sigma, rule=2)$y
      probe.scale <- 1 - scaling.x*(1 - theoretic.d.s.all) # equals  scaling.x*theoretic.d.s.all + (1-scaling.x)

      ## It is mandatory to limit the scaling factor - scalings <0 confound correction (like hook)
      probe.scale <- pmax(0.1, probe.scale)
      probe.scale <- pmin(1.0, probe.scale)

      ## Update pm and mm data
      pmOrMm(affyData)[,chip.idx] <- pmOrMm(affyData)[,chip.idx] / probe.scale
      
      if (plot.images) {
        means <- getGeneralDegradationMeans(x[subset.probes], as.vector(pmOrMm(affyData[,chip.idx]))[subset.probes],
                                            min.x, max.x)
        PlotDx(means, min.x, max.x, scale.x)
        tongs.cor <- ComputeTongs(x, as.vector(pmOrMm(affyData[,chip.idx])), sigma, probeset.id, min.x, max.x)
        PlotTongs(tongs.cor)
      }
    }

    stats.all[chip.idx, ] <- c(unlist(c(theoretic.d.s,hook.params)), decay)
  }
  
  return(new("AffyDegradationBatch", location.type = location.type,
             afbatch = affyData, stats = stats.all,
             means.pm = means.pm.all, means.mm = means.mm.all))
}

## Settings for index/absolue correction
kCorrectionTypeSettings <- list(
   index    = list(type.x = "index",    scale.x = 1,  min.x = 1, max.x = 11,  n.breaks = 11),
   absolute = list(type.x = "distance", scale.x = 60, min.x = 0, max.x = 600, n.breaks = 12))


setGeneric("d", function(object) standardGeneric("d"))
setMethod("d", "AffyDegradationBatch", function(object) 
          structure(object@stats[,"decay"], names = rownames(object@stats)))


setGeneric("afbatch", function(object) standardGeneric("afbatch"))
setMethod("afbatch", "AffyDegradationBatch", function(object) object@afbatch)


setGeneric("plotDx", function(object) standardGeneric("plotDx"))
setMethod("plotDx", "AffyDegradationBatch", function(object) {
  means = object@means.pm 
  x = as.numeric(rownames(means))

  with(kCorrectionTypeSettings[[ object@location.type ]], { ## Use index/absolute specific min.x, max.x, scale.c

    plot(0, xlim=c(min.x, max.x), ylim=c(0.1,1.2), pch=".", col="grey20",
         xlab="x", ylab="D(x)",
         panel.first = grid())

    cols <- rainbow(ncol(means))
    for (i in seq(1, ncol(means))) {
      d.x <- means[, i] / means[1, i]
      theoretic.d <- TryToFitDecayFunction((x - min.x)/scale.x, d.x, d.x[length(d.x)], 0.6) # - min.x to scale to 0
  
      points(x, d.x, pch=17, col=cols[i])
      xx <- seq(min(0, min.x), max.x, length.out = 100)
      lines(xx, DecayFunction((xx - min.x) / scale.x, theoretic.d$d1, theoretic.d$d2, theoretic.d$shiftY),
            lty = 2, lwd = 2, col = cols[i]) 
    }
    legend("topright", colnames(means), fill=cols)
  })
})


## Computes average intensity per interval or probe index
getGeneralDegradationMeans <- function(x, intensity, min.x = 1, max.x = 11, breaks = 11) {
  subset.x <- x >= min.x & x <= max.x
  x <- x[subset.x]
  intensity <-  intensity[subset.x]

  if(length(table(x)) <= breaks) { # small number of positions -> index correction -> simple average
    intensity <- log10(intensity)
    mean.intensity <- tapply(intensity, x, mean)
    return(data.frame(x = as.numeric(names(mean.intensity)), i = 10^(mean.intensity)))  
  } else {
    whist <- hist(x, breaks = breaks, plot=F)
    intervals<-length(whist$counts)

    absoluteMeans <- data.frame(x = rep(0,intervals), i = rep(0,intervals))
    for(i in 1:intervals){  
      ## if (whist$counts[i] > 1) {
      absoluteMeans[i,"x"] <- mean(c(whist$breaks[i], whist$breaks[i+1]))
      absoluteMeans[i,"i"] <- 10^mean(log10(intensity[x >= whist$breaks[i] & x < whist$breaks[i+1]]), na.rm = T)
    }
    return(absoluteMeans)
  }
}


# Function to plot decay function d(x) for x in k,L
PlotDx <-  function(means, min.x = 1, max.x = 11, scale.x = 1, ...) {
  d.x <- means$i / means$i[1]
  theoretic.d <- TryToFitDecayFunction((means$x - min.x)/scale.x, d.x, d.x[length(d.x)], 0.6) # - min.x to scale to 0
  
  plot(0, xlim=c(0,max.x), ylim=c(0.1,1.2), pch=".", col="grey20",
       xlab="x", ylab="D(x)",
       panel.first = grid(), ...) #           main=paste("Absolute probe location bias for",currentSample))
  points(means$x, d.x, pch=17)

  xx <- seq(min(0, min.x), max.x, length.out = 100)
  lines(xx, DecayFunction((xx - min.x) / scale.x, theoretic.d$d1, theoretic.d$d2, theoretic.d$shiftY),
        lty = 2, col = "grey60", lwd = 2)
}


ProbesetLogAverage <- function(x) { mean(log10(x)) }
ProbesetAverage <- function(x) { mean((x)) }
## ProbesetLogAverage <- function(x) { tukey.biweight(log10(x)) }
## ProbesetAverage <- function(x) { tukey.biweight((x)) }
