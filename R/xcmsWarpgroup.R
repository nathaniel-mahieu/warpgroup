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
#' @param sc.aligned.lim Integer. Peak bounds after alignment are considered the same if they are within this limit.
#' @param output.groups Boolean. If \code{TRUE} the output is a list of warpgroup outputs for every group rather than an xcmsSet.  Allows for better integration parameter selection with \code{\link{warpgroupsToXs}}.
#' @param sc.aligned.factor Float. Experimental feature where graph edges are weighted proportionally to the distance between the aligned peak bounds. Higher numbers emphasize closer peak bounds.
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
  sc.aligned.lim,
  output.groups=F,
  sc.aligned.factor = 1,
  detailed.groupinfo = F,
  min.peaks = 1
) { #Handles parallelization
  cat("This is a wrapper for the warpgroup algorithm to make it compatible with XCMS. The warpgroup algorithm performs peak grouping/clustering between samples, finds consensus peak bounds which describe similar regions of a peak for each group, and finds those consensus bounds in samples where a peak was not detected. The default output of this wrapper is an xcmsSet object for compatibility, but a list of groups, each containing the warpgroup determined peak bounds can be obtained by setting output.groups=T.\n")
  
  xs = add.raw.sc(xs)
  
  tryCatch(redisIncrBy("countTotal", length(xs@groupidx)), error=function(e) NULL)
  
  groups = foreach(
    params = iter.gwparams(xs, xr.l, rt.max.drift, ppm.max.drift),
    .packages = c("warpgroup", "dtw", "igraph"),
    .errorhandling = "pass",
    .noexport = c("xr.l"),
    .inorder=T
  ) %dopar% {
    tryCatch(redisDecrBy("countTotal", 1), error=function(e) NULL)
    
    groups = tryCatch(
{
  
  groups = warpgroup(
    ps = params$ps, # sc, scmin, scmax, sample
    eic.mat.s = params$eic.mat.s,
    
    sc.max.drift,
    sc.aligned.lim,
    sc.aligned.factor,
    detailed.groupinfo,
    min.peaks,
    tw="dtw"
  )
  
  groups = llply(groups, function(x) {
    
    rts = laply(seq(nrow(x)), function(i) {
      params$eic.mat[x[i,"sample"],,"rt"][floor(x[i,c("sc", "scmin", "scmax")])]
      })
    colnames(rts) = c("rt", "rtmin", "rtmax")
    
    cbind(x, rts, pn = params$gidx[x[,"n"]])
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


eic.mat.from.list = function(eic.l, length.target = Inf, upsample.force = F) {
  length = max(sapply(eic.l, function(x) {length(x$intensity)}))
  length.target.n = length.target
  
  if (length.target > 0 & length.target < 1) { #cat("Downsampling longest EIC to", length.target * 100, "percent.  This is", length*length.target, "scans."); 
                                               length.target.n = length*length.target }
  if (length.target > 1 & length.target < Inf) { #cat("Downsampling longest EIC to", length.target, "scans.  This is", round(length.target/length * 100), "percent of the longest EIC."); 
                                                 length.target.n = length.target }
  if (length.target.n > length & upsample.force == F) { #cat("Length target would upsample the timeseries. Instead, adjusting length target to longest timeseries:", length, "override with upsample.force=T."); 
                                                        length.target.n = length }
  
  
  eic.mat.s = array(integer(1), dim=c(length(eic.l), length.target.n, 2), dimnames=list(NULL, NULL, c("rt", "intensity")))
  for (i in seq(eic.l)) {
    eic.inter = approx(eic.l[[i]]$rt, eic.l[[i]]$intensity, n = length.target.n)
    eic.mat.s[i,,] = do.call(cbind, eic.inter)
    }
  
  eic.mat.s
}

iter.gwparams = function(xs, xr.l, rt.max.drift, ppm.max.drift) {
  it <- iter(xs@groupidx)
  
  nextEl = function() {
    gidx <- nextElem(it)  # throws "StopIteration"
    
    ps = xs@peaks[gidx,,drop=F]
    
    maxrt = min(sapply(xr.l, function(x) { max(x@scantime) }))
    
    rts = c(
      min(ps[,"rtmin"])-rt.max.drift, 
      max(ps[,"rtmax"])+rt.max.drift
    )
    rts[rts <= 0] = 1
    rts[rts > maxrt] = maxrt
    start.rt = min(rts)   
    
    scans = laply(xr.l, function(x) {
      h = sapply(rts, function(y) {which.min(abs(x@scantime - y))})
      c(h, diff(h))
      })
    
    #Per sample, summed EIC matrix for alignment
    mzmin = min(ps[,"mzmin"] - ps[,"mzmin"] * ppm.max.drift/1E6)
    mzmax = max(ps[,"mzmax"] + ps[,"mzmax"] * ppm.max.drift/1E6)
    
    eic.l = lapply(seq(xr.l), function(i) {
      l = rawEIC(xr.l[[i]], mzrange = c(mzmin, mzmax), scanrange = scans[i,c(1,2)])
      l$rt = xr.l[[i]]@scantime[l$scan]
      l
      })
    
    eic.mat = eic.mat.from.list(eic.l)
    eic.mat.s = aperm(eic.mat[,,2])
    
    ps.m = ps[,c("sc", "scmin", "scmax", "sample"),drop=F]
    for (r in seq(nrow(ps))) {
      s = ps[[r,"sample"]]
      for (i in seq(ps[r,c("sc", "scmin", "scmax")])) {
        ps.m[r,c("sc", "scmin", "scmax")][i] = which.min(abs(ps[r,c("rt", "rtmin", "rtmax")][i] - eic.mat[s,,"rt"]))
      }
    }
    
    list(
      eic.mat.s = eic.mat.s,
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