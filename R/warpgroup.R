warpgroup = function(
  ps, # sc, scmin, scmax, sample
  eic.mat.s,
  sc.max.drift,
  sc.aligned.lim
) {
  ps = cbind(ps, n=1:nrow(ps))
  
  #Find warped alignment between sample summed EICs
  dtw.step = buildStepMat(eic.mat.s) #dtw.mat = buildDtwMat(eic.mat.s)
  
  
  #Find predicted scan end and begin in all warp spaces.
  sc.warps = array(numeric(), dim = c(nrow(ps), nrow(ps), 3), dimnames=list(NULL,NULL, c("sc", "scmin", "scmax")))
  sc.warps.diffs = sc.warps
  for(i in seq_along(ps[,1])) { 
    a = ps[i,];
    for(j in seq_along(ps[,1])) { 
      b= ps[j,];
      a.sc.a = dtw.step[[a["sample"]]][[b["sample"]]](
        a[c("sc", "scmin", "scmax")]
      )
      b.sc.a = b[c("sc", "scmin", "scmax")]
      sc.warps[i,j,] = a.sc.a
      sc.warps.diffs[i,j,] = (a.sc.a - b.sc.a)
    }
  }
  
  
  #Define Grouping based on alignments
  if (nrow(ps) < 2) {
    cm.mem = c("1" = 1)
  } else {
    match.mat = aaply(abs(sc.warps.diffs), c(1,2), 
                      function(sc.ds) (
                        sc.ds[2] < sc.aligned.lim & sc.ds[3] < sc.aligned.lim |
                          sc.ds[1] < sc.aligned.lim*.75 & sc.ds[2] < sc.aligned.lim |
                          sc.ds[1] < sc.aligned.lim*.75 & sc.ds[3] < sc.aligned.lim
                      )
    )
    topology = which(match.mat, arr.ind=T)
    g2 = graph.data.frame(topology, directed=F)
    cm = walktrap.community(g2)
    cm.mem = membership(cm)
    
    g.l = split(cm.mem, cm.mem)
  }
  
  
  #Define peak bounds!
  cb.l = llply(unique(cm.mem), function(x) {
    pns.g = as.numeric(names(cm.mem[cm.mem == x]))
    sc.d = sc.warps.diffs[pns.g,pns.g,,drop=F] # How many scans apart the centwave peak bounds are after DTW
    sc.a = sc.warps[pns.g,pns.g,,drop=F]
    
    if(length(pns.g) > 1) {
      zscores = aperm(
        aaply(sc.d, c(3), scale),
        c(2,3,1)
      )
      
      voters = zscores > -1 & zscores < 1
      sc.a[!voters] = NA
      
      p.sc.params = aaply(sc.a, c(3), colMedians, na.rm=T)
    } else {
      p.sc.params = aaply(sc.a, c(3), .drop=F, identity)
    }
    
    colnames(p.sc.params) = pns.g
    p.sc.params
  })
  
  
  #Fillpeaks: find peakbounds pf missing peaks with with consensus-bounds found above and warps between raw data
  mb.l = llply(cb.l, function(cb) {
    found = ps[as.numeric(colnames(cb)), "sample"]
    all.samp = seq(ncol(eic.mat.s))
    missing = which(!all.samp %in% found)
    
    if (length(missing) < 1) return(NULL)
    
    mb = aperm(
      aaply(missing, 1, .drop=F, function(j) {
        aaply(colnames(cb), 1, .drop=F, function(i) {
          dtw.step[[ps[as.numeric(i),"sample"]]][[j]](
            cb[,i]
          )
        })
      })
    )
    dimnames(mb) = c(dimnames(cb), list(missing))
    
    foo = aaply(mb, c(1,3), median, .drop=F)
    ns = dimnames(foo)[1:2]
    dim(foo) = dim(foo)[-3]
    dimnames(foo) = ns
    foo
  })
  
  
  # Return Data
  llply(seq(cb.l), function(i) {
    cb = cb.l[[i]]
    cb.r = rbind(cb, n = as.numeric(colnames(cb)), sample = ps[as.numeric(colnames(cb)),"sample"])
    peaks = cb.r
    
    mb = mb.l[[i]]
    if (!is.null(mb)) {
      mb.r = rbind(mb, sample = as.numeric(colnames(mb)))
      
      peaks = merge(cb.r, mb.r, "row.names", all=T)
      rownames(peaks) = peaks[,"Row.names"]
      peaks = as.matrix(peaks[,-1])
    }
    
    peaks = aperm(peaks)[,c("n", "sample", "sc", "scmin", "scmax")]
    rownames(peaks) = NULL
    
    peaks[order(peaks[,"sample"]),]
  })
}