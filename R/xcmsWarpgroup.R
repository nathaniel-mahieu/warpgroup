group.warpgroup = function(
  xs,
  xr.l,
  sc.max.drift, 
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
    params = iter.gwparams(xs, xr.l, sc.max.drift, ppm.max.drift),
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
    min.peaks
  )
  
  groups = llply(groups, function(x) {
    x[,c("sc", "scmin", "scmax")] = x[,c("sc", "scmin", "scmax")] + params$start.scan
    
    cbind(x, pn = params$gidx[x[,"n"]])
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


iter.gwparams = function(xs, xr.l, sc.max.drift, ppm.max.drift) {
  it <- iter(xs@groupidx)
  
  nextEl = function() {
    gidx <- nextElem(it)  # throws "StopIteration"
    
    ps = xs@peaks[gidx,,drop=F]
    
    maxscan = min(sapply(xr.l, function(x) { length(x@scantime) }))
    
    scans = c(
      min(ps[,"scmin"])-sc.max.drift, 
      max(ps[,"scmax"])+sc.max.drift
    )
    scans[scans <= 0] = 1
    scans[scans > maxscan] = maxscan
    start.scan = min(scans)-1
    
    ps.m = ps[,c("sc", "scmin", "scmax", "sample"),drop=F]; ps.m[,c("sc", "scmin", "scmax")] = ps.m[,c("sc", "scmin", "scmax")] - start.scan;
    
    
    #Per sample, summed EIC matrix for alignment
    mzmin = min(ps[,"mzmin"] - ps[,"mzmin"] * ppm.max.drift/1E6)
    mzmax = max(ps[,"mzmax"] + ps[,"mzmax"] * ppm.max.drift/1E6)
    
    eic.mat.s = matrix(integer(1), ncol=length(xr.l), nrow=length(scans[1]:scans[2]))
    for (i in seq_along(xr.l)) { 
      eic = rawEIC(xr.l[[i]], mzrange = c(mzmin, mzmax), scanrange = scans)$intensity
      eic.mat.s[1:length(eic),i] = eic
    }
    
    list(
      eic.mat.s = eic.mat.s,
      ps = ps.m,
      gidx = gidx,
      start.scan = start.scan
    )
  }
  
  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
}



add.raw.sc = function(xs) {
  xs@peaks = xs@peaks[,-which(colnames(xs@peaks) %in% c("sc", "scmin", "scmax")),drop=F]
  
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