plot_peaks_bounds = function(eic, bounds) {
  eic.mat.s = apply(eic, 2, function(x) { max = max(x); if (max == 0) max = 1; x/max })
  
  df = melt(eic.mat.s)
  colnames(df) = c("Scan", "Sample", "Intensity")
  
  bounds = cbind(bounds, sample.peak = seq_along(bounds[,1]))
  
  int = array(logical(), dim=dim(eic.mat.s)); int[] = F
  for (r in seq(nrow(bounds))) int[floor(bounds[r,"scmin"]):floor(bounds[r,"scmax"]), bounds[r,"sample"]] = T
  df.int = data.frame(df, Integrated = melt(int)[,3])
  
  df.bounds = melt(data.frame(bounds), id.vars = c("sample.peak", "sample"), measure.vars=c("scmin", "scmax"))
  colnames(df.bounds) = c("sample.peak", "Sample", "Bound", "Scan")
  Intensity = aaply(as.matrix(df.bounds[,c("Sample", "Scan")]), 1, function(x) {
    eic.mat.s[x["Scan"],x["Sample"]]
  })
  df.bounds = cbind(df.bounds, Intensity)
  
  ggplot(df.int, aes(y = Intensity, x = Scan, group=factor(Sample))) + 
    geom_line(data = df.int, colour="#a0a0a0") +
    geom_line(data=subset(df.int, Integrated==T), mapping=aes(group=factor(Sample)), colour="#D63647", size=1) + 
    geom_point(data=df.bounds, colour="#D63647", size=4) +
    geom_text(data = df.bounds[!duplicated(df.bounds[,"sample.peak"]),], mapping = aes(label = Sample), x = max(df.int$Scan)+5) +
    theme(
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(colour="#F2EBEC")
    ) + facet_wrap(~Sample, scales="free_y")
}


