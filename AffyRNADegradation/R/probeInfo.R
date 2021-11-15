## Collects probe index information and returns data table in the same order as AffyBatch intensity matrix
getProbeInfo.index <- function(affyData) {
  probeset.id <- probeNames(affyData)

  ## For each set, probe are indexed in reverse order of their
  ## occurence in the intensity matrix
  ## exception: different ordering for PrimeView CDFs
  if (cdfName(affyData) != "PrimeView") {
    probe.index <- ave(seq(length(probeset.id), 1, -1), probeset.id, FUN=rank)
  } else {
    probe.index <- ave(seq(1, length(probeset.id),  1), probeset.id, FUN=rank)
  }

  return(data.frame(probeset.id = probeset.id, x = probe.index, stringsAsFactors=F))
}

## Helper function to return appopriate probe info
getProbeInfo <- function(affyData, location.type = "index", location.file.dir = NULL) {
  if (location.type == "index") {
    return(getProbeInfo.index(affyData))
  } else {
    return(getProbeInfo.absolute(affyData, location.file.dir))
  }
}

## Collect absolute probe location information
## and returns data table in the same order as AffyBatch intensity matrix
getProbeInfo.absolute <- function(affyData, location.file.dir = NULL) {
  location.file = file.path(location.file.dir, paste(cdfName(affyData), ".probe_distances.Rd", sep=""))
  
  tryCatch({ suppressWarnings(load(location.file)) },
           error = function(err) {
             stop(paste("Cannot open probe location file. Make sure that directory given by",
                        " location.file.dir (\"", location.file.dir, "\") contains the respective file",
                        " for the current chip type (", cdfName(affyData), "). Probe location files",
                        " for all Affymetrix expression arrays can be found at",
                        " http://www.izbi.uni-leipzig.de/downloads_links/programs/rna_integrity.php", sep=""),
                  call.=FALSE)
           })

  ## Add affy indices as rownames to that table (xy2indices)
  options(scipen = 99) # This is important: otherwise indices like 400000 are 
                       # transformed into "4e5" leading to bugs 
  rownames(probeDists) <- apply(probeDists[c("Probe.X","Probe.Y")], 1,
            function(x) xy2indices(x[1], x[2], abatch=affyData))

  ## Create table in same order as pm() containing probeset IDs, distance 
  ## ...  mapping the probe distances using the affy probe numbers (x*width +y)
  ## @note Control probes must be handled carefully: insert pseudo-values
  probeInfo <- data.frame(pmI = pm(affyData[,1]), probeset.id = "control",
    x = NA, stringsAsFactors=F)
  probeInfo[,1] <- NULL # delete row with intensities - only used as a template
  probeInfo[rownames(probeDists), c("probeset.id","x")] <- probeDists[,c("Probe.Set.Name","Probe.Distance")]

  return(probeInfo)
}
  
