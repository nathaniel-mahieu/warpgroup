warpgroup = function(
  ps, # mz, sc, scmin, scmax, sample
  eic.mat.s,
  start.scan,
  
  sc.drift.lim,
  ppm.lim,
  sc.lim
) {
  ps = cbind(ps, n=1:nrow(ps))
  
  #Find warped alignment between sample summed EICs
  dtw.step = buildStepMat(eic.mat.s)
  
  
  #Find predicted scan end and begin in all warp spaces.
  sc.warps = array(numeric(), dim = c(nrow(ps), nrow(ps), 4), dimnames=list(NULL,NULL, c("sc", "sc.min", "sc.max", "sc.wmean")))
  sc.warps.diffs = sc.warps
  for(i in seq_along(ps[,1])) { 
    a = ps[i,];
    for(j in seq_along(ps[,1])) { 
      b= ps[j,];
      a.sc.a = dtw.step[[a["sample"]]][[b["sample"]]](
        a[c("sc", "sc.min", "sc.max", "sc.wmean")] - start.scan
      )
      b.sc.a = b[c("sc", "sc.min", "sc.max", "sc.wmean")] - start.scan
      sc.warps[i,j,] = a.sc.a
      sc.warps.diffs[i,j,] = (a.sc.a - b.sc.a)
    }
  }
  
  
  #Define Grouping based on alignments
  if (nrow(ps) < 2) {
    cm.mem = c("1" = 1)
  } else {
    #ppm.mat = outer(ps[,"mz"], ps[,"mz"], function(x,y) abs(x-y)/x*1E6)
    #ppm.mat[upper.tri(ppm.mat, T)] = NA
    #topo.ppm = melt(ppm.mat)
    #topo.ppm = topo.ppm[!is.na(topo.ppm[,"value"]),]
    #topo.ppm[,"value"] = topo.ppm[,"value"] * 2
    #topo.ppm[topo.ppm[,"value"] < 1,"value"] = 1
    
    #g2.ppm = graph.data.frame(topo.ppm, directed=F)
    #cm.ppm = walktrap.community(g2.ppm, weights=1/topo.ppm$value)
    #plot(cm.ppm, g2.ppm, edge.width = 1/topo.ppm$value)
    
    match.mat = aaply(abs(sc.warps.diffs), c(1,2), 
                      function(sc.ds) (sc.ds[1] < sc.lim | sc.ds[4] < sc.lim) | (sc.ds[2] < sc.lim & sc.ds[3] < sc.lim)
    )
    topology = which(match.mat, arr.ind=T)
    g2 = graph.data.frame(topology, directed=F)
    cm = walktrap.community(g2)
    cm.mem = membership(cm)
  }
  
  
  #Define peak bounds!
  sc.params.l = llply(unique(cm.mem), function(x) {
    pns.g = as.numeric(names(cm.mem[cm.mem == x]))
    sc.d = sc.warps.diffs[pns.g,pns.g,,drop=F]
    sc.a = sc.warps[pns.g,pns.g,,drop=F]
    
    if(length(pns.g) > 1) {
      pct90 = aaply(abs(sc.a), c(3), quantile, 0.9) # Should this be sc.aad?
      
      voters = laply(seq(dim(sc.a)[3]), function(i) sc.a[,,i,drop=F] <= pct90[[i]])
      voters = aperm(voters, c(2,3,1))
      
      sc.a.vote = sc.a; sc.a.vote[!voters] = NA
      p.sc.params = aaply(sc.a.vote, c(3), colMeans, na.rm=T)
    } else {
      p.sc.params = aaply(sc.a, 3, .drop=F, function(x) x)
    }
    colnames(p.sc.params) = pns.g
    p.sc.params
  })
  p.sc.params = do.call("cbind", sc.params.l)
  consensus.bounds = aperm(p.sc.params)
  
  
  #Fillpeaks: find missing peakbounds with with consensus peak bounds and warps
  all.samp = seq(xr.l)
  pgs = merge(ps, data.frame(cm.mem, n=names(cm.mem)), by = "n", all=T)
  new.gs = split(cm.mem, cm.mem)
  
  miss.samps = llply(new.gs, function(x) which(!all.samp %in% subset(pgs, cm.mem == x[1])[,"sample"]))
  
  fill.samps = llply(seq(miss.samps), function(i) {
    miss = miss.samps[[i]]  
    pres = names(new.gs[[i]])
    consensus.bounds[pres,]
    fills = array(numeric(), dim=c(length(miss), length(pres), ncol(consensus.bounds)))
    dimnames(fills) = list(miss,pres,colnames(consensus.bounds))
    
    for(i in miss) { 
      for(j in pres) {          
        fills[as.character(i),j,] = dtw.step[[pgs[j,"sample"]]][[i]](consensus.bounds[j,])
      }
    }
    foo = aaply(fills, c(1,3), mean)
    dim(foo) = dim(fills)[c(1,3)]
    colnames(foo) = dimnames(fills)[[3]]
    rownames(foo) = dimnames(fills)[[1]]
    foo
  })
  
  
  #Return data
  llply(seq(new.gs), function(i) {
    pnames = names(new.gs[[i]])
    found = cbind(
      n = as.numeric(pnames), 
      consensus.bounds[pnames, c("sc", "sc.min", "sc.max"),drop=F] + start.scan, 
      samp = ps[as.numeric(pnames),"sample"],
      modularity = NA,
      warp.cost = NA
    )
    
    n.miss = nrow(fill.samps[[i]])
    missed = cbind(
      n = rep(NA, n.miss), 
      fill.samps[[i]][, c("sc", "sc.min", "sc.max"),drop=F] + start.scan, 
      samp = as.numeric(rownames(fill.samps[[i]])),
      modularity = rep(NA, n.miss),
      warp.cost = rep(NA, n.miss)
    )
    
    foo = rbind(found,missed)
  })
}