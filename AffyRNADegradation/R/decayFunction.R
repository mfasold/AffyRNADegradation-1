DecayFunction <- function(x, d.inf, lambda.inv, shiftY) { (1 - d.inf)*exp(-lambda.inv*x) + d.inf + shiftY }

TryToFitDecayFunction <- function(x, means, defaultD1, defaultD2, defaultShiftY= 0) {
  params <- list(d1 = defaultD1, d2 = defaultD2, shiftY = defaultShiftY)
  
  p = try({
    p2 = try({ fit = nls(means ~ DecayFunction(x, d1, d2, shiftY),
                 start=list(d1 = defaultD1, d2 = defaultD2, shiftY = defaultShiftY),
                 lower = c(0, 0.01, -0.2), upper = c(1.1, 5, 1), algorithm="port")}, silent=T) #Use lower and upper bounds for robustness

    ## If the fit does not converge, assume d1=0 
    if(inherits(p2, "try-error")) {
      cat("Warning: Convergence error. Setting d1=0.\n")
      fit <- nls(means ~ DecayFunction(x, 0, d2, shiftY),
        start = list(d2 = defaultD2, shiftY = defaultShiftY),
        lower = c(0.005, -0.2), upper = c(2, 0.3), algorithm="port")      
      params$d1 <- 0
    } else {
      params$d1 <- coef(fit)["d1"]   
    }
    
    params$shiftY <- coef(fit)["shiftY"]
    params$d2 <- coef(fit)["d2"]
  }, silent=T)
  if(inherits(p, "try-error")) {
    cat("Convergence error. Using default degradation parameters. \n")
  }
  return(params)
}

