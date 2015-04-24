plotGroup.xs = function(i, xs, xr.l, sc.pad = 20, mz.pad = 0) {
  g = xs@groupidx[[i]]
  ps = xs@peaks[g,,drop=F]
  
  scan.range = c(floor(min(ps[,"scmin"], na.rm=T))-sc.pad, ceiling(max(ps[,"scmax"], na.rm=T))+sc.pad)
  eic.mat = matrix(numeric(), ncol = nrow(ps), nrow = scan.range[2]-scan.range[1]+1)
  
  for (j in seq(ps[,1])) {
    mz.range = unname(c(ps[j,"mzmin"]-mz.pad, ps[j,"mzmax"]+mz.pad))
    eic.mat[,j] = rawEIC(xr.l[[ps[j,"sample"]]], mzrange = mz.range, scanrange = scan.range)$intensity 
  }
  
  plot.me = melt(eic.mat)
  colnames(plot.me) = c("scan", "peak", "int")
  plot.me[,"peak"] = factor(plot.me[,"peak"])
  
  annos = data.frame(ps)
  annos[,"n"] = seq(annos[,1])
  colnames(annos)[c(9, 11, 12, 13)] = c("peak","scwm", "scmax", "scmin")
  annos[,c(11, 12, 13)] = annos[,c(11, 12, 13)] - min(scan.range)
  
  ggplot() + 
    geom_path(data = plot.me, aes(y=int, x= scan, colour=peak)) + 
    geom_point(data=annos, mapping=aes(y=0, x=scmax), colour="#000000") + 
    geom_point(data=annos, mapping=aes(y=0, x=scmin), colour="#000000") +
    geom_vline(data=annos, mapping=aes(xintercept=scwm), colour="#000000", alpha=0.3)+
    facet_grid(peak ~ ., scales="free_y")
}
