dtwFunc = function(v1, v2, keep=F) {
  dtw(
    v1,
    v2,
    keep=keep, 
    open.end = T,
    open.begin = T, 
    step.pattern = asymmetricP05, 
    window.type = "none"
  )
}

prepEicMat = function(eic.mat, n.pad=0) {
  pad.mat = matrix(0, ncol = ncol(eic.mat), nrow = n.pad)
  
  eic.mat = filter(eic.mat, c(1,1,1)/3)
  eic.mat[is.na(eic.mat)] = 0 

  eic.mat = aaply(eic.mat, 1, function(x) { 
    m = max(x); m[m==0] = 1;
    x/m
    })
  
  eic.mat = rbind(pad.mat, eic.mat, pad.mat)
  eic.mat
}

buildDtwMat = function(eic.mat) {
  n.pad = floor(nrow(eic.mat) * .25)
  eic.mat = prepEicMat(eic.mat, n.pad)
  
  llply(seq_along(eic.mat[1,]), function(i) {
    llply(seq_along(eic.mat[1,]), function(j) {
      dtwFunc(eic.mat[,i], eic.mat[,j], keep=T)
    })
  })
}

buildStepMat = function(eic.mat) {  
  n.pad = floor(nrow(eic.mat) * .25)
  eic.mat = prepEicMat(eic.mat, n.pad)
  
  llply(seq_along(eic.mat[1,]), function(i) {
    llply(seq_along(eic.mat[1,]), function(j) {
      align = dtwFunc(eic.mat[,i], eic.mat[,j])
      
      stepfun(
        align$index1[-1] - n.pad,
        align$index2 - n.pad
      )
    })
  })
}

merge.peaks = function(x) {
  c = count(x[,"samp"])
  dups = c[c[,"freq"] > 1, "x"]
  
  if(length(dups) < 1) return(x)
  
  merge = llply(dups, function(y) which(x[,"samp"] == y))
  
  replacements = laply(merge, function(y) {
    x2 = x[y,]
    
    cbind(
      sc = mean(x2[,"sc"]), 
      sc.min = min(x2[,"sc.min"]),
      sc.max = max(x2[,"sc.max"]),
      centWave.pn = 0,
      mzmin = min(x2[,"mzmin"]),
      mzmax = max(x2[,"mzmax"]),
      mz = mean(x2[,"mz"]),
      samp = x2[1,"samp"],
      modularity = x2[1,"modularity"],
      warp.cost = x2[1,"warp.cost"]
    )
    
  })
  rbind(x[-unlist(merge),], replacements)
} 

add.raw.sc = function(xs) {
  rt.scans = matrix(NA, ncol=3, nrow=nrow(xs@peaks), dimnames=list(NULL, c("sc", "sc.min", "sc.max")))
  
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