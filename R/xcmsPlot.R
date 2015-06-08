plotGroup.xs = function(i, xs, xr.l, sc.pad = 20, mz.pad = 0, type=1) {
  g = xs@groupidx[[i]]
  ps = xs@peaks[g,,drop=F]
  
  ps = ps[order(ps[,"sample"]),]
  
  eic.mat = makeEicMat(ps, xr.l, sc.pad = sc.pad, mz.pad = mz.pad)
  scan.range = range(c(ps[,"scmin"]))
  
  bounds = ps[,c("scmin","scmax", "sample")]; bounds[,1:2] = bounds[,1:2] - min(scan.range) +1 + sc.pad
  bounds[,"sample"] = seq(nrow(bounds))
  
  plot_peaks_bounds(eic.mat, na.omit(bounds), type=type)
}

makeEicMat = function(ps, xr.l, sc.pad = 20, mz.pad = 0) {
  #sc.start = min(ps[,"scmin"], na.rm=T) + 1
  
  maxscan = max(sapply(xr.l, function(x) { length(x@scantime) }))
  
  scan.range = c(floor(min(ps[,"scmin"], na.rm=T))-sc.pad, ceiling(max(ps[,"scmax"], na.rm=T))+sc.pad)
  scan.range[scan.range < 1] = 1; scan.range[scan.range > maxscan] = maxscan
  eic.mat = matrix(numeric(), ncol = nrow(ps), nrow = scan.range[2]-scan.range[1]+1)
  
  mz.rangeg = unname(c(min(ps[,"mzmin"]-mz.pad, na.rm=T), max(ps[,"mzmax"]+mz.pad, na.rm=T)))
  
  for (j in seq(ps[,1])) {
    if (is.na(ps[j,"mzmin"])) {mz.range = mz.rangeg} else {
    mz.range = unname(c(ps[j,"mzmin"]-mz.pad, ps[j,"mzmax"]+mz.pad))
    }
    eic.mat[,j] = rawEIC(xr.l[[ps[j,"sample"]]], mzrange = mz.range, scanrange = scan.range)$intensity 
  }
  
  eic.mat
}