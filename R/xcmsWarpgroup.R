#' A wrapper for supplying xcmsSet objects to Warpgroup.
#' 
#' \code{group.warpgroup} returns an xcmsSet with the original groups and peaks replaced by those returned by the \code{\link{warpgroup}} algorithm.
#' 
#' This function simply iterates over each group in the xcmsSet and runs \code{\link{warpgroup}} on them. The resulting warpgroup output is then re-integrated and returned as a new xcmsSet.
#' The resulting xcmsSet is ready for a call to \code{\link{diffreport}} and does not need to be grouped.  The results will be redundant and should be filtered appropriately.
#' 
#' @param xs An xcmsSet object with initial grouping information.
#' @param xr.l A list contanining xcmsRaw objects for each sample in the xcmsSet object in order.
#' @param sc.max.drift Integer.  The maximum time drift expected for a peak in the data set in scans.  Used when looking for missing peaks.
#' @param ppm.max.drift Integer. The maximum mass drift expected for a peak in the data set in ppm.  Used when looking for missing peaks.
#' @param rt.aligned.lim Integer. Peak bounds after alignment are considered the same if they are within this limit.
#' @param eic.resample.target If less than one the resulting EICs will be of length max*length.target.  If greater than 1 resulting EICs will be of length length.target.  If Inf resulting EICs will be of length max.
#' @param smooth.n Number of points to consider for the moving average smoothing of each EIC. 1 disables smoothing.
#' @param normalize Boolean.  If TRUE all EICs will be normalized to 1.
#' @param output.groups Boolean. If \code{TRUE} the output is a list of warpgroup outputs for every group rather than an xcmsSet.  Allows for better integration parameter selection with \code{\link{warpgroupsToXs}}.
#' @param sc.aligned.factor Float. Experimental feature where graph edges are weighted proportionally to the distance be   tween the aligned peak bounds. Higher numbers emphasize closer peak bounds.
#' @param detailed.groupinfo Boolean. Returns several extra descriptors of the warping and graph clustering.
#' @param min.peaks Integer. Groups with fewer peaks than \code{min.peaks} are skipped.
#' 
#' @return An xcmsSet with groups and peaks replaced by warpgroup generated ones.
#' 
#' @seealso See \url{https://github.com/nathaniel-mahieu/warpgroup} for examples.

group.warpgroup = function(
  xs,
  xr.l,
  rt.max.drift, 
  ppm.max.drift, 
  rt.aligned.lim,
  smooth.n,
  normalize = T,
  eic.resample.target = Inf,
  output.groups=F,
  sc.aligned.factor = 1,
  detailed.groupinfo = F,
  min.peaks = 1
) { #Handles parallelization
  cat("This is a wrapper for the warpgroup algorithm to make it compatible with XCMS. The warpgroup algorithm performs peak grouping/clustering between samples, finds consensus peak bounds which describe similar regions of a peak for each group, and finds those consensus bounds in samples where a peak was not detected.\nIf peak detection was poor it is likely that many peak groups will be similar or redundant.\n\nWarpgroup operates on the assumption that all samples will have similar EIC traces. As such warpgroup will incorrectly define peak groups in cases where samples are very different.")
  
  xs = add.raw.sc(xs)
  
  tryCatch(redisIncrBy("countTotal", length(xs@groupidx)), error=function(e) NULL)
  
  groups = foreach(
    params = iter.gwparams(xs, xr.l, rt.max.drift, ppm.max.drift, eic.resample.target, smooth.n, normalize),
    .packages = c("warpgroup", "dtw", "igraph"),
    .errorhandling = "pass",
    .noexport = c("xr.l"),
    .inorder=T
  ) %dopar% {
    tryCatch(redisDecrBy("countTotal", 1), error=function(e) NULL)
    
    groups = tryCatch(
{
  #params = nextElem(params)
  
  sc.aligned.lim = round(rt.aligned.lim /  mean(diff(params$eic.mat[1,,"rt"]),na.rm=T))
  
  groups = warpgroup(
    ps = params$ps, # sc, scmin, scmax, sample
    eic.mat.s = aperm(params$eic.mat[,,"intensity"]),

    sc.aligned.lim = sc.aligned.lim,
    sc.aligned.factor = sc.aligned.factor,
    detailed.groupinfo = detailed.groupinfo,
    min.peaks = min.peaks,
    tw="dtw"
  )
  
  groups = llply(groups, function(x) {
    x[, c("sc", "scmin", "scmax")][x[,c("sc", "scmin", "scmax")] < 1] = 1
    x[, c("sc", "scmin", "scmax")][x[,c("sc", "scmin", "scmax")] > length(params$eic.mat[1,,1])] = length(params$eic.mat[1,,1])
    
    rts = laply(seq(nrow(x)), function(i) {
      params$eic.mat[x[[i,"sample"]],,"rt"][floor(x[i,c("sc", "scmin", "scmax")])]
      })
    colnames(rts) = c("rt.raw", "rtmin.raw", "rtmax.raw")
    
    cbind(x[,-which(colnames(x) %in% c("sc", "scmin", "scmax"))], rts, pn = params$gidx[x[,"n"]])
  })
  
},
error=paste,
finally = return
    )

groups
  }

if (output.groups) return(groups)
return(warpgroupsToXs(xs, groups, xr.l)) # Merge back into traditional xcms workflow
}


iter.gwparams = function(xs, xr.l, rt.max.drift, ppm.max.drift, eic.resample.target, smooth.n, normalize) {
  it <- iter(xs@groupidx)
  maxrt = min(sapply(xr.l, function(x) { max(x@scantime) }))
  
  nextEl = function() {
    gidx <- nextElem(it)  # throws "StopIteration"
    ps = xs@peaks[gidx,,drop=F]
    
    # Define RT region for EIC generation
    rts = c(
      min(ps[,"rtmin"])-rt.max.drift, 
      max(ps[,"rtmax"])+rt.max.drift
    )
    rts[rts <= 0] = 1
    rts[rts > maxrt] = maxrt

    scans = laply(xr.l, function(x) {
      sapply(rts, function(y) {which.min(abs(x@scantime - y))})
      })
    
    # Define m/z region for EIC generation
    mzmin = min(ps[,"mzmin"] - ps[,"mzmin"] * ppm.max.drift/1E6)
    mzmax = max(ps[,"mzmax"] + ps[,"mzmax"] * ppm.max.drift/1E6)
    
    eic.l = lapply(seq(xr.l), function(i) {
      l = rawEIC(xr.l[[i]], mzrange = c(mzmin, mzmax), scanrange = scans[i,c(1,2)])
      l$rt = xr.l[[i]]@scantime[l$scan]
      do.call(rbind, l)
      })
    
    eic.mat = eicMatFromList(eic.l, eic.resample.target = eic.resample.target, smooth.n = smooth.n, normalize = normalize)
    
    ps.m = ps[,c("sc", "scmin", "scmax", "sample"),drop=F]
    for (r in seq(nrow(ps))) {
      s = ps[[r,"sample"]]
      for (i in seq(ps[r,c("sc", "scmin", "scmax")])) {
        ps.m[r,c("sc", "scmin", "scmax")][i] = which.min(abs(ps[r,c("rt", "rtmin", "rtmax")][i] - eic.mat[s,,"rt"]))
      }
    }
    
    list(
      eic.mat = eic.mat,
      ps = ps.m,
      gidx = gidx
    )
  }
  
  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
}



add.raw.sc = function(xs) {
  xs@peaks = xs@peaks[,!colnames(xs@peaks) %in% c("sc", "scmin", "scmax"),drop=F]
  
  rt.scans = matrix(NA, ncol=3, nrow=nrow(xs@peaks), dimnames=list(NULL, c("sc", "scmin", "scmax")))
  
  for (i in unique(xs@peaks[,"sample"])) {
    which.ps = which(xs@peaks[,"sample"]==i)
    rt.sc = stepfun(
      xs@rt$corrected[[i]][-1],
      seq_along(xs@rt$corrected[[i]])
    )
    rt.scans[which.ps,] = rt.sc(xs@peaks[which.ps,c("rt", "rtmin", "rtmax")])
  }
  xs@peaks = cbind(xs@peaks, rt.scans)
  xs
}