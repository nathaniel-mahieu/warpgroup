warpgroupsToXs = function(xs, groups, xr.l, ppm.padding=1) {
  cat("Converting warpgroups to xcmsSet.\nNote: The xcmsSet returned by this function does not need fillpeaks().\nCaution: diffreport() performs further processing on the peak groups before reporting statistics. Specifically it discards overlapping groups. This could remove groups which describe different portions of a peak found by the warpgrouping yet overlap.  If this behavior is not desired statstics can easily be performed on the raw warpgroup data retrieved by setting output.groups=T.\n")
  
  group.l = unlist(groups, F)
  bad.gs = unique(c(
    which(!sapply(group.l, is.matrix))
  ))
  if (length(bad.gs) > 0) warning(paste(length(bad.gs), "groups were removed due to errors."))
  for (bg in rev(bad.gs)) group.l[[bg]] = NULL
  
  pt.l = foreach(
    params = iter.integrateparams(group.l, xs, xr.l, ppm.padding),
    .errorhandling = "pass",
    .noexport = c("xr.l"),
    .inorder=T
  ) %do% {
    ints = integrate.simple(
      scanmat.l = params$scanmat,
      eic.l = params$eic.l
    )
    
    cbind(params$g, ints)
  }
  
  pt.l = lapply(seq(pt.l), function(i) cbind(pt.l[[i]], new.gidx = i))
  pt = do.call("rbind", pt.l)
  xs@peaks = pt[order(pt[,"sample"]),]
  xs = buildGroups(xs, xs@peaks[,"new.gidx"])
  
  return(xs)
}

iter.integrateparams = function(group.l, xs, xr.l, ppm.padding) {
  it <- iter(group.l)
  
  nextEl = function() {
    g <- nextElem(it)  # throws "StopIteration"
    
    g = g[!duplicated(g[,"sample"]),,drop=F]
    
    x = xs@peaks[g[,"pn"],,drop=F]
    
    mzrange.g = c(
      min(x[,"mzmin"], na.rm=T) - max(x[,"mzmin"], na.rm=T) * ppm.padding / 1E6, 
      max(x[,"mzmax"], na.rm=T) + max(x[,"mzmax"], na.rm=T) * ppm.padding / 1E6
    )
    
    foreach(i=seq(nrow(g))) %do% {
      p = x[i,]
      pg = g[i,]
      
      mzrange = mzrange.g
      if (!is.na(p["mzmin"])) mzrange = unname(c(p["mzmin"] - p["mzmin"] * ppm.padding / 1E6, p["mzmax"] + p["mzmin"] * ppm.padding / 1E6))
      
      scanrange = as.numeric(c(
        floor(pg["scmin"]), 
        ceiling(pg["scmax"])
      ))
      scanrange[scanrange < 1] = 1
      maxscan = length(xr.l[[pg["sample"]]]@scantime)
      scanrange[scanrange > maxscan] = maxscan 
      eic = rawEIC(xr.l[[pg["sample"]]], mzrange = mzrange, scanrange = scanrange)      
      eic = cbind(eic, rt=xs@rt$corrected[eic$scan])
      
      scans = foreach(i=scanrange[1]:scanrange[2]) %do% getScan(xr.l[[pg["sample"]]], i, mzrange)
      scanmat = do.call("rbind", scans)
      
      list(
        eic = eic,
        scanmat = scanmat,
        g = g
      )
    }
  }
  
  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
}

integrate.simple = function(scanmat.l, eic.l) {
  int.mat = matrix(numeric(), nrow=length(scanmat.l), ncol=8, dimnames=list(NULL, c("mz", "rt", "rt.half", "into", "maxo", "mzmin", "mzmax", "rtmin", "rtmax")))
  
  foreach(i=seq(scanmat.l)) %do% {
    scan.mat = scanmat.l[[i]]
    eic = eic.l[[i]]
    
    min = quantile(scan.mat[,"intensity"], .5)
    mz = mean(scan.mat[scan.mat[,"intensity"] > c(min), "mz"])
    
    weightedi = diff(eic$rt) * sapply(
      2:length(eic$intensity), function(x){
        mean(eic$intensity[c(x - 1, x)])
      })
    into = sum(weightedi)
    rtwmean = sum(weightedi * eic$rt[-1]) / into
    
    rt.half = mean(eic$rt[1], eic$rt[length(eic$rt)])
    
    int.mat[i,] = c(
      mz,
      rtwmean,
      rt.half,
      into,
      max(eic$intensity),
      mzrange[[1]],
      mzrange[[2]],
      min(rts),
      max(rts)
    )
  }
  
  int.mat
}

buildGroups = function(xs, pgs) {
  if(nrow(xs@peaks) != length(pgs)) stop("Peaks number does not equal group assignments.")
  
  #Make groupidx
  pgs = as.numeric(pgs)
  groups = llply(unique(pgs), .progress="text", function(x) which(x == pgs))
  xs@groupidx = groups
  
  #Fill groupval
  classnames = as.character(unique(xs@phenoData$class))
  groupval = matrix(numeric(), nrow=length(xs@groupidx), ncol=7+length(classnames), dimnames=list(NULL, c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax","npeaks",classnames)))
  
  for(i in seq_along(xs@groupidx)) {
    g = xs@groupidx[[i]]
    ps = xs@peaks[g,,drop=F]
    groupval[i,] = 
      c(
        median(ps[,"mz"], na.rm=T),
        range(ps[,"mz"], na.rm=T),
        median(ps[,"rt"], na.rm=T),
        range(ps[,"rt"], na.rm=T),
        length(g),
        rep(NA, length(classnames))
      )
  }
  xs@groups = groupval
  xs
}
