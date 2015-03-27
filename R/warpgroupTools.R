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