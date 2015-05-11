plot.warpgroup = function(group, eic.mat.s, scan.start) {
  ps = group 
  
  long.plot = melt(eic.mat.s)
  colnames(long.plot) = c("scan", "sample", "intensity")
  
  colors = llply(ps[,"sample"], function(n) cbind(sample=n, scan = floor(ps[ps[,"sample"]==n,"sc.min"]):ceiling(ps[ps[,"sample"]==n,"sc.max"]) - start.scan, color=n))
  colors = do.call("rbind", colors)
  
  lp2 = merge(long.plot, colors, by=c("sample", "scan"), all=T)
  
  ggplot() + 
    geom_path(data = lp2, mapping = aes(x=scan, y= intensity), colour="#a3a3a3") + 
    geom_path(data = subset(lp2, !is.na(color)), mapping = aes(x=scan, y= intensity, colour = factor(sample)))   + 
    facet_wrap(~sample, scales="free_y") + theme_bw()+ theme(legend.position="none")
}


plot.warpgroup = function(group, xs, xr.l, sc.pad = 20, mz.pad = 0.001) {
  g = group
  ps = xs@peaks[g[,"pn"],,drop=F]
  ps[,c("sc","scmin","scmax","sample")] = g[,c("sc","scmin","scmax","sample")]
  
  eic.mat = makeEicMat(ps, xr.l, sc.pad, mz.pad)

  g[,c("sc", "scmin", "scmax")] = g[,c("sc", "scmin", "scmax")] - min(g[,"scmin"]) + 1 + sc.pad
  plot_peaks_bounds(eic.mat, g)
}


plot_peaks_bounds = function(eic.mat.s, bounds) {
  eic.mat.s = apply(eic.mat.s, 2, function(x) { x/max(x) })
  eic.mat.s = sapply(seq(ncol(eic.mat.s)), function(i) { eic.mat.s[,i]+1.1*(i-1) })
  
  df = melt(eic.mat.s)
  colnames(df) = c("Scan", "Sample", "Intensity")
  
  int = array(logical(), dim=dim(eic.mat.s)); int[] = F
  for (r in seq(nrow(bounds))) int[floor(bounds[r,"scmin"]):floor(bounds[r,"scmax"]), bounds[r,"sample"]] = T
  df.int = melt(int)
  colnames(df.int) = c("Scan", "Sample", "Integrated")
  
  df.int = data.frame(df, Integrated = melt(int)[,3])
  
  df.bounds = melt(data.frame(bounds), id.vars = "sample", measure.vars=c("scmin", "scmax"))
  colnames(df.bounds) = c("Sample", "Bound", "Scan")
  Intensity = aaply(as.matrix(df.bounds[,c("Sample", "Scan")]), 1, function(x) {
    eic.mat.s[x["Scan"],x["Sample"]]
    })
  df.bounds = cbind(df.bounds, Intensity)
  
  if(F) {
  ggplot(df.int, aes(y = Intensity, x = Scan)) + 
    geom_line() +
    geom_line(data=subset(df.int, Integrated==T), mapping=aes(colour="red")) +
    geom_point(data=df.bounds, aes(y=1)) + 
    facet_wrap(~Sample, ncol=2, scales="free_y") +
    theme(legend.position="none", axis.text.y=element_blank(), axis.ticks=element_blank()) +
    theme(#axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
}
  
  ggplot(df.int, aes(y = Intensity, x = Scan, group=factor(Sample))) + 
    geom_line(colour="#a0a0a0") +
    geom_line(data=subset(df.int, Integrated==T), mapping=aes(group=factor(Sample)), colour="#D63647", size=1) + 
    geom_point(data=df.bounds, colour="#D63647", size=4) +
    theme(
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour="#F2EBEC")
      )
}

plot.warpgroup = function(prepath=".", ps, cm.mem, eic.mat, eic.mat.s, p.sc.params) {
  
  foo = merge(ps, data.frame(cm.mem, n=names(cm.mem)), by = "n", all=T)    
  
  colnames(eic.mat) = ps[,"n"]
  long.plot = melt(eic.mat)
  colnames(long.plot) = c("scan", "n", "intensity")
  
  pbs = cbind(t(as.matrix(p.sc.params)), n=as.numeric(colnames(p.sc.params)))
  colnames(pbs) = c("sc", "scmin", "scmax", "scwm", "n")
  pbs = merge(pbs, ps[,c("sample", "n")], by = "n")
  
  colors = llply(pbs[,"n"], function(n) cbind(n=n, scan = floor(pbs[pbs[,"n"]==n,"scmin"]):ceiling(pbs[pbs[,"n"]==n,"scmax"]), color=n))
  colors = do.call("rbind", colors)
  
  lp2 = merge(long.plot, foo[,c("n", "cm.mem")], by="n", all=T)
  lp3 = merge(lp2, colors, by=c("n", "scan"), all=T)
  lp4 = merge(lp3, pbs[,c("n","sample")], by = "n", all=T)
  
  ggplot() + 
    geom_path(data = lp4, mapping = aes(x=scan, y= intensity), colour="#a3a3a3") + 
    geom_path(data = subset(lp4, !is.na(color)), mapping = aes(x=scan, y= intensity, colour = factor(cm.mem))) + 
    geom_vline(data=data.frame(pbs), mapping=aes(xintercept=scwm), colour="#000000", alpha=0.1)+
    geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmin), colour="#000000")+
    geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmax), colour="#000000")+
    facet_wrap(~n) + xlim(100,350)+ theme_bw()+ theme(legend.position="none")
  
  ggplot() + 
    geom_path(data = lp4, mapping = aes(x=scan, y= intensity), colour="#a3a3a3") + 
    geom_path(data = subset(lp4, !is.na(color)), mapping = aes(x=scan, y= intensity), colour="blue") + 
    #geom_vline(data=data.frame(pbs), mapping=aes(xintercept=scwm), colour="#000000", alpha=0.1)+
    #geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmin), colour="#000000")+
    #geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmax), colour="#000000")+
    facet_wrap(~n)+ theme_bw()+ theme(legend.position="none")
  
  
  
  ggplot(data.frame(ps), aes(x=sc, y=mz, colour=factor(sample))) + geom_point() + geom_segment(mapping=aes(y=mz, yend=mz, x=sc.min, xend=sc.max)) #+ geom_segment(mapping=aes(y=mzmin, yend=mzmax, x=rt, xend=rt))
  #ggplot(data.frame(cbind(ps, aligned.scans, n = 1:nrow(ps))), aes(x=sc.a, y=mz, colour=factor(sample), label=n)) + geom_point()+ geom_segment(mapping=aes(y=mz, yend=mz, x=sc.min.a, xend=sc.max.a)) + geom_segment(mapping=aes(y=mzmin, yend=mzmax, x=sc.a, xend=sc.a)) + geom_text(mapping=aes(y=mz+.002, x= sc.a+0.2))
  ggplot(melt(eic.mat.s), aes(y=value, x= Var1, colour=factor(Var2))) + geom_path() + facet_wrap(~Var2)
  ggplot(melt(eic.mat), aes(y=value, x= Var1, colour=factor(Var2))) + geom_path() + facet_grid(Var2 ~ ., scales="free_y") 
  
  
  
  dev.off()
  #plot(aligns[[1]][[2]],type="twoway")
  #plot(aligns[[2]],type="threeway")
  #plotRawMzRtI(ps, xr.l)
  
  long.plot = melt(eic.mat)
  colnames(long.plot) = c("scan", "n", "intensity")
  
  pbs = cbind(t(as.matrix(p.sc.params)), n=as.numeric(colnames(p.sc.params)))
  colnames(pbs) = c("sc", "scmin", "scmax", "scwm", "n")
  pbs = merge(pbs, ps[,c("sample", "n")], by = "n")
  
  colors.orig = llply(ps[,"n"], function(n) cbind(n=n, scan = floor(ps[ps[,"n"]==n,"sc.min"]):ceiling(ps[ps[,"n"]==n,"sc.max"])-start.scan, color=n))
  colors.orig = do.call("rbind", colors.orig)
  
  lp2.o = merge(long.plot, foo[,c("n", "cm.mem")], by="n", all=T)
  lp3.o = merge(lp2.o, colors.orig, by=c("n", "scan"), all=T)
  lp4.o = merge(lp3.o, ps[,c("n","sample")], by = "n", all=T)
  
  ggplot() + 
    geom_path(data = lp4.o, mapping = aes(x=scan, y= intensity), colour="#a3a3a3") + 
    geom_path(data = subset(lp4.o, !is.na(color)), mapping = aes(x=scan, y= intensity), colour="blue") + 
    #geom_vline(data=data.frame(pbs), mapping=aes(xintercept=scwm), colour="#000000", alpha=0.1)+
    #geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmin), colour="#000000")+
    #geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmax), colour="#000000")+
    facet_wrap(~n) + xlim(100,350)+ theme_bw()+ theme(legend.position="none")
  
  #Peaks colored by grouping old peak bounds
  ggplot() + 
    geom_path(data = lp4.o, mapping = aes(x=scan, y= intensity), colour="#a3a3a3") + 
    geom_path(data = subset(lp4.o, !is.na(color)), mapping = aes(x=scan, y= intensity, colour = factor(cm.mem))) + 
    #geom_vline(data=data.frame(pbs), mapping=aes(xintercept=scwm), colour="#000000", alpha=0.1)+
    #geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmin), colour="#000000")+
    #geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmax), colour="#000000")+
    facet_wrap(~n) + xlim(100,350)+ theme_bw()+ theme(legend.position="none")
  
}