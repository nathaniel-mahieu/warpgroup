
if(F){
  #Plots
  foo = merge(ps, data.frame(cm.mem, n=names(cm.mem)), by = "n", all=T)    
  
  eic.mat.plot = laply(seq_along(eic.mat[1,]), function(i) {
    x = eic.mat[,i]
    y = x
    x[] = NA
    x[foo[i,"sc.min"]:foo[i,"sc.max"] - min(scans)+1] = y[foo[i,"sc.min"]:foo[i,"sc.max"] - min(scans)+1]
    x
  })
  long.plot = melt(eic.mat.plot)
  long.plot.2 = merge(long.plot, foo, by.x="Var1", by.y="n", all.x=T)
  long.plot.2 = long.plot.2[order(long.plot.2[,"Var2"]),]
  
  pbs = cbind(t(as.matrix(p.sc.params)), Var1=as.numeric(colnames(p.sc.params)))
  colnames(pbs) = c("sc", "scmin", "scmax", "scwm", "Var1") 
  ggplot() + 
    geom_path(data = melt(t(eic.mat)), mapping = aes(x=Var2, y= as.numeric(value)), colour="#ffffff") + 
    geom_path(data = long.plot.2, aes(x = Var2, y = as.numeric(value), colour = factor(cm.mem))) + 
    geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmin), colour="#000000")+
    geom_point(data=data.frame(pbs), mapping=aes(y=0, x=scmax), colour="#000000")+
    geom_vline(data=data.frame(pbs), mapping=aes(xintercept=scwm), colour="#000000", alpha=0.1)+
    facet_grid(Var1~.)
  
  
  
  
  #ggplot(data.frame(cbind(ps, aligned.scans,pn = 1:nrow(ps))), aes(x=rt, y=mz, colour=factor(sample))) + geom_point() + geom_segment(mapping=aes(y=mz, yend=mz, x=rtmin, xend=rtmax)) + geom_segment(mapping=aes(y=mzmin, yend=mzmax, x=rt, xend=rt))
  #ggplot(data.frame(ps), aes(x=rt, y=mz, colour=factor(sample))) + geom_point() + geom_segment(mapping=aes(y=mz, yend=mz, x=rtmin, xend=rtmax)) + geom_segment(mapping=aes(y=mzmin, yend=mzmax, x=rt, xend=rt))
  #ggplot(data.frame(cbind(ps, aligned.scans, n = 1:nrow(ps))), aes(x=sc.a, y=mz, colour=factor(sample), label=n)) + geom_point()+ geom_segment(mapping=aes(y=mz, yend=mz, x=sc.min.a, xend=sc.max.a)) + geom_segment(mapping=aes(y=mzmin, yend=mzmax, x=sc.a, xend=sc.a)) + geom_text(mapping=aes(y=mz+.002, x= sc.a+0.2))
  #dtwPlot(aligns[[2]][[4]],type="density")
  #plot(aligns[[2]][[4]], type="threeway")
  ggplot(melt(eic.mat.s), aes(y=value, x= Var1, colour=factor(Var2))) + geom_path() + facet_grid(Var2 ~ .)
  ggplot(melt(eic.mat), aes(y=value, x= Var1, colour=factor(Var2))) + geom_path() + facet_grid(Var2 ~ ., scales="free_y") 
  
  
  #plot(aligns[[1]][[2]],type="twoway")
  #plot(aligns[[2]],type="threeway")
  #plotRawMzRtI(ps, xr.l)
  
  
  #Add psg information
  psgs = subset(data.frame(xs@peaks), psg %in% ps[,"psg"])
  
  mz.rnd = round(psgs[,"mz"], 2)
  mz.rnd.u = unique(mz.rnd)
  psg.u = unique(ps[,"psg"])
  psg.presence = matrix(F, ncol = length(psg.u), nrow = length(mz.rnd.u), dimnames=list(as.character(mz.rnd.u),as.character(psg.u)))
  for (i in seq_along(mz.rnd)) {
    psg.presence[as.character(mz.rnd[[i]]), as.character(psgs[i,"psg"])] = T
  }
  
  hist(rowSums(psg.presence))
  psg.presence[rowSums(psg.presence)>4,]
  
  
  p.sc.a.mat = array(numeric(), dim = c(nrow(psgs),nrow(psgs), 4))
  for(i in seq_along(ps[,1])) { 
    a = ps[i,]
    for(j in seq_along(ps[,1])) { 
      b= ps[j,]
      align = aligns[[a["sample"]]][[b["sample"]]]
      sc.sc = stepfun(
        align$index1[-1],
        align$index2
      )
      
      p.sc.a.mat[i,j,] = isGroupMatch(a,b,sc.sc, sc.lim = 5, scans = scans)
    }
  }

}


buildSinglePsg = function(xs) {
  psgs = llply(split(xs, seq_along(xs@phenoData$class)), function(xs.1) {
    xsa.1 = xsAnnotate(xs.1)
    xsa.1 = groupFWHM(xsa.1)
    xsa.1 = groupCorr(xsa.1)
    
    psgs = melt(xsa.1@pspectra)
    psgs[order(psgs[,"value"]), "L1"]
  })
  psgs.m = melt(psgs)
  psg.f = as.numeric(factor(paste(psgs.m[,1], psgs.m[,2], sep=".")))
  psg.f
}
