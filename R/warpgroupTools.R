dtwFunc = function(v1, v2, keep=T) {
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

getStep = function(a, tw, n.pad) {
  if (tw == "dtw") {
    return(stepfun(a$index1[-1] - n.pad, a$index2 - n.pad))
  } else if (tw == "ptw") {
    return(stepfun(a$index1[-1] - n.pad, a$index2 - n.pad))
  }
}

buildTwMat = function(eic.mat, pct.pad = 0, tw="dtw") {
  n.pad = floor(nrow(eic.mat) * pct.pad)
  eic.mat = prepEicMat(eic.mat, n.pad)
  n = ncol(eic.mat)
  cor = list()
  cor$tw = tw
  cor$twmat = llply(seq(n), function(x) {vector("list", n) })
  cor$stepmat = llply(seq(n), function(x) {vector("list", n) })
  cor$npad = n.pad
  
  if (tw == "dtw") {
    twFunc = dtwFunc
  } else if (tw == "ptw") {
    twFunc = ptwFunc
  }
  
  for (i in seq(n)) {
    for (j in seq(n)) {
      cor$twmat[[i]][[j]] = twFunc(eic.mat[,i], eic.mat[,j])
      cor$stepmat[[i]][[j]] = getStep(cor$twmat[[i]][[j]], cor$tw, cor$npad)
    }
  }

  cor
}