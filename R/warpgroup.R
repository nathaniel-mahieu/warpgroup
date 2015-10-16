#' Warpgroup a set of peak bounds.
#' 
#' \code{warpgroup} returns a list of all distinct chromatographic regions determined from the supplied peak bounds for every sample.
#' 
#' This function implements the time-warping based chromatogram correspondance determination, splitting of supplied peak bounds into distinct groups, and selection of the correct region to integrate for missing values.
#' The warpgroup approach can be applied to a variety of situations but is most effective when supplied groups are of high quality.
#' 
#' Warpgroup first determines the correspondance between supplied EICs using a time warping technique (by default dynamic time warping).  Based on this correspondance warpgroup splits the supplied peak bounds into sub groups, each which describe similar chromatographic regions.  Warpgroup then determines the appropriate integration region for samples missing peaks.
#' 
#' @param ps Matrix. Rows representing peak regions detected in \code{eic.mat.s} and containing columns \code{sc, scmin, scmax, sample}.
#' @param eic.mat.s Matrix. Chromatograms, one per sample. Columns corresponding to samples and rows corresponding to time with values of intensity.
#' @param sc.aligned.lim Integer. Peak bounds after alignment are considered the same if they are within this limit.
#' @param detailed.groupinfo Logical. Returns several extra descriptors of the warping and graph clustering.
#' @param sc.aligned.factor Float. Experimental feature where graph edges are weighted proportionally to the distance between the aligned peak bounds. Higher numbers emphasize closer peak bounds.
#' @param min.peaks Integer. Supplied groups with fewer than this number of peak bounds will be skipped.
#' @param tw Character. \code{"tw"} or \code{"dtw"} Allows the selection of the desired correspondance algorithm
#' @param pct.pad Float. Defines how much padding is added to the ends of each chromatogram.
#' 
#' @return A list, one entry for each distinct chromatographic region detected.  Within each entry is a matrix, with one row for each sample describing the start and end of the corresponding peak region as indices of the supplied \code{eic.mat.s}.
#' 
#' @seealso See \url{https://github.com/nathaniel-mahieu/warpgroup} for examples.

warpgroup = function(
  ps, # sc, scmin, scmax, sample
  eic.mat.s,
  sc.aligned.lim,
  sc.aligned.factor = 1,
  detailed.groupinfo = F,
  min.peaks = 1,
  tw = "dtw",
  pct.pad = 0.1
) {
  if (nrow(ps) < min.peaks) stop("Less peaks in starting group than min.peaks.")
  if (nrow(eic.mat.s) < 3) stop("EIC contains too few points.")
  
  ps = cbind(ps, n=1:nrow(ps))
  
  #Find warped alignment between sample summed EICs
  tw.l = buildTwList(eic.mat.s, tw = tw, pct.pad = pct.pad)
  
  #Find predicted scan end and begin in all warp spaces.
  sc.warps = array(numeric(), dim = c(nrow(ps), nrow(ps), 3), dimnames=list(NULL,NULL, c("sc", "scmin", "scmax")))
  sc.warps.diffs = sc.warps
  for(i in seq_along(ps[,1])) { 
    a = ps[i,];
    for(j in seq_along(ps[,1])) { 
      b= ps[j,];
      a.sc.a = tw.l[[a[["sample"]]]][[b[["sample"]]]]$step(
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
    g.characteristics.l = list(matrix(NA, ncol=5, nrow=1, dimnames = list(NULL, c("group.degree", "group.density", "group.eccentricity", "group.closeness", "group.parent.modularity"))))
  } else {
    match.mat = aaply(abs(sc.warps.diffs), c(1,2), function(sc.ds) (
      sc.ds[2] <= sc.aligned.lim & sc.ds[3] <= sc.aligned.lim
      )
    )

    match.wts = aaply(abs(sc.warps.diffs), c(1,2), 
      function(sc.ds) 1/(1.1^((sc.ds[2] + sc.ds[3])/sc.aligned.factor))
      )

    topology = which(match.mat, arr.ind=T); top.wts = topology
    top.wts = cbind(topology, weight = match.wts[topology])

    g = graph.data.frame(top.wts, directed=F)
    cm = walktrap.community(g)
    cm.mem = membership(cm)
    g.l = split(cm.mem, cm.mem)
    
    g.characteristics.l = llply(unique(cm.mem), function(g) {
      g = g.l[[as.character(g)]]
      
      topo2 = subset(data.frame(top.wts, row.names=NULL), 
                     X1 %in% names(g) & X2 %in% names(g))
      g2 = graph.data.frame(topo2, directed=F)
      
      cbind(
        group.degree = round(degree(g2)/(vcount(g2)*2),2),
        group.density = graph.density(g2),
        group.eccentricity = eccentricity(g2),
        group.closeness = closeness(g2),
        group.parent.modularity = modularity(cm, g)
      )
    })
  }
  
  
  #Define peak bounds!
  cb.l = llply(unique(cm.mem), function(x) {
    pns.g = as.numeric(names(cm.mem[cm.mem == x]))
    sc.d = sc.warps.diffs[pns.g,pns.g,,drop=F] # How many scans apart the centwave peak bounds are after DTW
    sc.a = sc.warps[pns.g,pns.g,,drop=F]
    
    if(length(pns.g) > 1) {
      # Pick voting peaks
      zscores = aperm(
        aaply(sc.d, c(3), function(x) {
          matrix(scale(c(x), center=F), ncol=ncol(x), byrow=F)
          }),
        c(2,3,1)
      )
      foo = aaply(zscores, c(1,3), function(x) {mean(abs(x))})
      
      rowtf = aaply(foo, 2, function(x) { 
        if (all(is.na(x))) {
          !is.na(x)
        } else if (min(x) < 1) {
          x <= 1
        } else { #Fallback if no peaks within 1 SD
          x < quantile(x, 0.75)
        }
      })
      rowtf2 = aperm(outer(rowtf, rep(T, ncol(rowtf)), "&"),c(2,3,1))
      celltf = abs(zscores) < 2
      
      #Iterate to converge chosen bounds
      for (y in 1:4) {
        sc.a[!rowtf2 | !celltf] = NA
        
        p.sc.params = round(aaply(sc.a, c(3), colMeans, na.rm=T))
      
        #Check peak bounds - should all be the same when projected into each others time domain.
        sc.warps.projected = array(numeric(), dim = dim(sc.a), dimnames=list(NULL,NULL, c("sc", "scmin", "scmax")))
        for (i in seq(pns.g)) {
          for (j in seq(pns.g)) {
            sc.warps.projected[i, j, ] = tw.l[[ps[pns.g[i],"sample"]]][[ps[pns.g[j],"sample"]]]$step(p.sc.params[,as.character(i)])
          }
        }
        
        sc.a = sc.warps.projected
      }
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
          tw.l[[ps[as.numeric(i),"sample"]]][[j]]$step(
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
  groups = llply(seq(cb.l), function(i) {
    cb = cb.l[[i]]
    cb.r = rbind(cb, n = as.numeric(colnames(cb)), sample = ps[as.numeric(colnames(cb)),"sample"], aperm(g.characteristics.l[[i]]))    
    peaks = cb.r
    
    mb = mb.l[[i]]
    if (!is.null(mb)) {
      mb.r = rbind(mb, sample = as.numeric(colnames(mb)))
      
      peaks = merge(cb.r, mb.r, "row.names", all=T)
      rownames(peaks) = peaks[,"Row.names"]
      peaks = as.matrix(peaks[,-1])
    }
    
    desired.columns = c("n", "sample", "sc", "scmin", "scmax")
    if(detailed.groupinfo) desired.columns = c(desired.columns, c("group.degree", "group.density", "group.eccentricity", "group.closeness", "group.parent.modularity"))
    
    peaks = aperm(peaks)[,desired.columns]
    rownames(peaks) = NULL
    
    peaks[order(peaks[,"sample"]),]
  })
  
  if(detailed.groupinfo & tw == "dtw") {
    groups = llply(groups, function(g) {
      d.phi.cum = matrix(numeric(), ncol=nrow(g), nrow = nrow(g))
      d.cum = matrix(numeric(), ncol=nrow(g), nrow = nrow(g))
      
      warp.consistency = array(numeric(), dim = c(nrow(g), nrow(g), 3))

      for (i in seq(nrow(g))) {
        for (j in seq(nrow(g))) {
          p = g[i,]
          
          sc = c(floor(p["scmin"]), ceiling(p["scmax"]))
          
          tw.m =  tw.l[[p["sample"]]][[g[j,"sample"]]]
          tw.m.rev =  tw.l[[g[j,"sample"]]][[p["sample"]]]
          
          warp.consistency[i,j,] = tw.m.rev$step(tw.m$step(c(sc, mean(sc)))) - c(sc,mean(sc))

          sc = sc + tw.m$npad
          if (length(dim(tw.m$path)) ==2) {
            indices = which.min(abs(tw.m$path[,1] - sc[1])):which.min(abs(tw.m$path[,1] - sc[2]))
          } else {
            indices=1:2
          }
          
          d.phi.cum[i,j] = (tw.m$d.phi[tail(indices, n=1)] - tw.m$d.phi[1]) / length(indices)
          d.cum[i,j] = sum(tw.m$d[indices]) / length(indices)
        }
      }
      
      diag(d.phi.cum)= NA
      diag(d.cum)= NA
      
      bestpeak = which.min(colMeans(aaply(warp.consistency, 3, rowSds)))
      
      cbind(g, dtw.distortion = rowMeans(d.phi.cum, na.rm=T), warped.distance = rowMeans(d.cum, na.rm=T), warp.consistency = rowMeans(abs(warp.consistency[,,3])))
      })
  }
  
  groups[
    which(sapply(groups, function(x) { sum(!is.na(x[!duplicated(x[,"sample"]),"n"])) >= min.peaks }))
    ]
}