plotGroup.xs = function(i, xs, xr.l, sc.pad = 20, mz.pad = 0) {
  g = xs@groupidx[[i]]
  ps = xs@peaks[g,,drop=F]
  
  ps = ps[order(ps[,"sample"]),]
  
  eic.mat = makeEicMat(ps, xr.l, sc.pad = sc.pad, mz.pad = mz.pad)
  scan.range = range(c(ps[,"scmin.raw"]))
  
  bounds = ps[,c("scmin.raw","scmax.raw", "sample")]; bounds[,1:2] = bounds[,1:2] - min(scan.range) +1 + sc.pad
  bounds[,"sample"] = seq(nrow(bounds))
  colnames(bounds) = c("scmin","scmax", "sample")
  
  plot_peaks_bounds(eic.mat, na.omit(bounds))
}

plot.warpgroup = function(group, xs, xr.l, sc.pad = 20, mz.pad = 0.001, type=1) {
  g = group
  ps = xs@peaks[g[,"pn"],,drop=F]
  ps[,c("sc","scmin","scmax","sample")] = g[,c("sc","scmin","scmax","sample")]
  
  eic.mat = makeEicMat(ps, xr.l, sc.pad, mz.pad)
  
  g[,c("sc", "scmin", "scmax")] = g[,c("sc", "scmin", "scmax")] - min(g[,"scmin"]) + 1 + sc.pad
  plot_peaks_bounds(eic.mat, g, type=type)
}

makeEicMat = function(ps, xr.l, sc.pad = 20, mz.pad = 0) {
  
  maxscan = min(sapply(xr.l, function(x) { length(x@scantime) }))
  
  scan.range = c(floor(min(ps[,"scmin.raw"], na.rm=T))-sc.pad, ceiling(max(ps[,"scmax.raw"], na.rm=T))+sc.pad)
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