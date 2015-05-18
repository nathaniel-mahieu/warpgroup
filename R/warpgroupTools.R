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

ptwFunc = function(v1, v2, keep=F) {
  library(ptw)
  ptw(
    v1,
    v2,
    init.coef=c(0,1,0),
    smooth.param = 0,
    trwdth = 20
  )
}

prepEicMat = function(eic.mat, n.pad=0) {
  pad.mat = matrix(0, ncol = ncol(eic.mat), nrow = n.pad)
  
  eic.mat.ma = filter(eic.mat, c(1,1,1)/3, sides=2)
  eic.mat.ma[is.na(eic.mat.ma)] = 0 
  eic.mat[] = eic.mat.ma@.Data
  
  eic.mat = t(aaply(eic.mat, 2, function(x) { 
    m = max(x); m[m==0] = 1;
    x/m
  }))
  
  eic.mat = rbind(pad.mat, eic.mat, pad.mat)
  eic.mat
}

getCorrespondance = function(v1, v2, tw = "dtw", n.pad = 0) {
  if (tw == "dtw") {
    dtw = dtwFunc(v1, v2)
    return(stepfun(
      dtw$index1[-1] - n.pad,
      dtw$index2 - n.pad
    ))
  } else if (tw == "ptw") {
    ptw = ptwFunc(v1, v2)
    return(stepfun(
      seq(ptw$warp.fun)[-1],
      ptw$warp.fun
      ))
  } else {
    stop("Unavailable time warping function selected.")
  }
}

buildDtwMat = function(eic.mat, pct.pad = 0.1, tw="dtw") {
  n.pad = floor(nrow(eic.mat) * pct.pad)
  eic.mat = prepEicMat(eic.mat, n.pad)
  n = seq(ncol(eic.mat))
  
  llply(n, function(i) {
    llply(n, function(j) {
      dtwFunc(eic.mat[,i], eic.mat[,j])
    })
  })
}

buildStepMat = function(eic.mat, pct.pad = 0.1, tw="dtw") {  
  n.pad = floor(nrow(eic.mat) * pct.pad)
  eic.mat = prepEicMat(eic.mat, n.pad)
  n = seq(ncol(eic.mat))
  
  llply(n, function(i) {
    llply(n, function(j) {
      getCorrespondance(eic.mat[,i], eic.mat[,j], tw="dtw", n.pad)
    })
  })
}