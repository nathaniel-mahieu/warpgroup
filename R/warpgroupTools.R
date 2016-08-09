#' A helper function to interpolate a list of EIC's into a rectangular array with optional down/up-sampling and smoothing
#' 
#' \code{eicMatFromList} returns an array containing each EIC's resulting intensity, retention time, and scan.
#' 
#' This function takes a list of EIC's (intensities, retention times, and indices for each scan) and interpolates them such that they are the same length and normalizes the intensities to 1.
#' 
#' @param eic.l A list of matrices respresnting EICs.  Each matrix has rows "rt", "intensity", and "scan"
#' @param eic.resample.target Numeric.  If less than one the resulting EICs will be of length max*length.target.  If greater than 1 resulting EICs will be of length length.target.  If Inf resulting EICs will be of length max.
#' @param upsample.force Boolean.  Must be true for resulting EICs to be of length > max.
#' @param smooth.n Integer. The number of points to include in the moving average.
#' @param normalize Boolean. If T all EICs will be normalized to 1.0
#' 
#' @return An array containing each EIC's resulting intensity, retention time, and scan.

eicMatFromList = function(eic.l, eic.resample.target = Inf, upsample.force = F, smooth.n = 1, normalize = T) {
  length = max(sapply(eic.l, function(x) {ncol(x)}))
  length.target.n = eic.resample.target
  
  if (eic.resample.target > 0 & eic.resample.target < 1) { #cat("Downsampling longest EIC to", eic.resample.target * 100, "percent.  This is", length*eic.resample.target, "scans."); 
    length.target.n = round(length*eic.resample.target) }
  if (eic.resample.target > 1 & eic.resample.target < Inf) { #cat("Downsampling longest EIC to", eic.resample.target, "scans.  This is", round(eic.resample.target/length * 100), "percent of the longest EIC."); 
    length.target.n = round(eic.resample.target) }
  if (length.target.n > length & upsample.force == F) { #cat("Length target would upsample the timeseries. Instead, adjusting length target to longest timeseries:", length, "override with upsample.force=T."); 
    length.target.n = length }
  
  
  eic.mat.s = array(integer(1), dim=c(length(eic.l), length.target.n, 3), dimnames=list(NULL, NULL, c("rt", "intensity", "scan")))
  for (i in seq(eic.l)) {
    eic = eic.l[[i]]
    interpolated.int = approx(eic["rt",], eic["intensity",], n = length.target.n)
    interpolated.scans = approx(eic["rt",], eic["scan",], n = length.target.n)
    
    x = filter(interpolated.int$y, rep(1, smooth.n)/smooth.n, sides=2)
    x[is.na(x)] = 0 
    x = x@.Data
    
    if (normalize)   m = max(x); m[m==0] = 1; x = x/m
    
    eic.mat.s[i,,c("rt")] = interpolated.int$x
    eic.mat.s[i,,c("intensity")] = x
    eic.mat.s[i,,c("scan")] = interpolated.scans$y
  }
  
  eic.mat.s
}



dtwFunc = function(v1, v2, keep=T) {
  dtw(
    v1,
    v2,
    keep=keep, 
    open.end = T,
    open.begin = T, 
    step.pattern = asymmetricP1, 
    window.type = "none"
  )
}

ptwFunc = function(v1, v2, keep=F) {
  library(ptw)
  ptw(
    v1,
    v2,
    init.coef=c(0,1,0),
    smooth.param = 0,
    trwdth = 20
  )
}

getStep = function(a, tw, n.pad) {
  if (tw == "dtw") {
    return(stepfun(a$index1[-1] - n.pad, a$index2 - n.pad))
  } else if (tw == "ptw") {
    return(stepfun(a$index1[-1] - n.pad, a$index2 - n.pad))
  }
}

buildTwList = function(eic.mat, pct.pad = 0, tw="dtw") {
  n.pad = floor(nrow(eic.mat) * pct.pad)
  
  pad.mat = matrix(0, ncol = ncol(eic.mat), nrow = n.pad)
  eic.mat = rbind(pad.mat, eic.mat, pad.mat)
  
  n = ncol(eic.mat)
  tw.l = rep(list(vector("list",n)), n)
  
  if (tw == "dtw") {
    twFunc <<- dtwFunc
  } else if (tw == "ptw") {
    twFunc <<- ptwFunc
  }
  
  for (i in seq(n)) {
    for (j in seq(n)) {
      #cor$twmat[[i]][[j]] = twFunc(eic.mat[,i], eic.mat[,j])
      warp = twFunc(eic.mat[,i], eic.mat[,j], keep=T)
      
      data = list()
      data$step = getStep(warp, tw, n.pad)
      data$path = cbind(warp$index1, warp$index2)
      data$d.phi = warp$costMatrix[ data$path ]
      data$d = warp$localCostMatrix[ data$path ]
      data$tw = tw
      data$npad = n.pad
      #data$warp = warp
      
      tw.l[[i]][[j]] = data
    }
    gc()
  }

  tw.l
}