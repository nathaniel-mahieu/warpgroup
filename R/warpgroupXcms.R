group.warpgroup = function(
  xs,
  xr.l=NULL,
  sc.max.drift, 
  ppm.max.drift, 
  sc.aligned.lim, 
  .parallel = F,
  output.groups=F
  ) { #Handles parallelization
  cat("This is a wrapper for the warpgroup algorithm to make it compatible with XCMS. The warpgroup algorithm performs peak grouping/clustering between samples, finds consensus peak bounds which describe similar regions of a peak for each group, and finds those consensus bounds in samples where a peak was not detected. The default output of this wrapper is an xcmsSet object for compatibility, but a list of groups, each containing the warpgroup determined peak bounds can be obtained by setting output.groups=T.\n")
  
  if (is.null(xr.l)) xr.l = llply(xs@filepaths, xcmsRaw, profstep=0)

  groups = llply(
    seq_along(xs@groupidx),
    .progress="text", 
    .parallel=.parallel, 
    .paropts=list(
      .packages=c("warpgroup", "dtw", "igraph", "xcms")#, .export = c("xs", "xr.l")
      ), 
    xs, xr.l,
    
    .fun = function(i, xs, xr.l) {
      cat("\r", i)
      tryCatch(
        xcms.warpgroup(
          xs = xs, 
          g.i = i, 
          xr.l = xr.l,
          sc.max.drift=sc.max.drift, 
          ppm.max.drift = ppm.max.drift, 
          sc.aligned.lim = sc.aligned.lim
          ),
        error=paste,
        finally = return
        )
      }
    )
  
  if (output.groups) return(groups)
  return(warpgroupsToXs(xs, groups, xr.l)) # Merge back into traditional xcms workflow
}

xcms.warpgroup = function(
  xs, 
  g.i, 
  xr.l,
  sc.max.drift, 
  ppm.max.drift, 
  sc.aligned.lim
  ) { #Handles xcms setup
  
  xs = add.raw.sc(xs)
  
  ps = xs@peaks[xs@groupidx[[g.i]],,drop=F] 
 
  scans = c(
    min(ps[,"scmin"])-sc.max.drift, 
    max(ps[,"scmax"])+sc.max.drift
    )
  scans[scans <= 0] = 1
  start.scan = min(scans)-1

  #Per sample, summed EIC matrix for alignment
  mzmin = min(ps[,"mzmin"] - ps[,"mzmin"] * ppm.max.drift/1E6)
  mzmax = max(ps[,"mzmax"] + ps[,"mzmax"] * ppm.max.drift/1E6)
  
  eic.mat.s = matrix(integer(1), ncol=length(xr.l), nrow=length(scans[1]:scans[2]))
  for (i in seq_along(xr.l)) { 
    eic = rawEIC(xr.l[[i]], mzrange = c(mzmin, mzmax), scanrange = scans)$intensity
    eic.mat.s[1:length(eic),i] = eic
  }
  
  ps.m = ps[,c("sc", "scmin", "scmax", "sample"),drop=F]; ps.m[,c("sc", "scmin", "scmax")] = ps.m[,c("sc", "scmin", "scmax")] - start.scan;
  
  groups = warpgroup(
    ps = ps.m, # sc, scmin, scmax, sample
    eic.mat.s = eic.mat.s,
    
    sc.max.drift,
    sc.aligned.lim
  )
  
  pn = as.matrix(data.frame(
    n = seq(xs@groupidx[[g.i]]), pn = xs@groupidx[[g.i]]
    ))
  llply(groups, function(g) {
    as.matrix(merge(g, pn, by="n", all.x = T))
    })
}




warpgroupsToXs = function(xs, groups, xr.l, ppm.padding=1) {
  cat("Converting warpgroups to xcmsSet.\nNote: The xcmsSet returned by this function does not need fillpeaks().\nCaution: diffreport() performs further processing on the peak groups before reporting statistics. Specifically it discards overlapping groups. This could remove groups which describe different portions of a peak found by the warpgrouping yet overlap.  If this behavior is not desired statstics can easily be performed on the raw warpgroup data retrieved by setting output.groups=T.\n")
  
  group.l = unlist(groups, F)
  bad.gs = which(!sapply(group.l, is.matrix))
  warning(paste(length(bad.gs), "groups were removed due to errors."))
  for (bg in rev(bad.gs)) group.l[[bg]] = NULL
  
  sum.dups = 0
  
  pt.l = llply(seq(group.l), .progress="text", function(i) {
    g = group.l[[i]]
    dupd = duplicated(g[,"sample"]); sum.dups = sum.dups + sum(dupd)
    g = g[!dupd,,drop=F]
    g = cbind(g, new.gidx = i)
    
    ints = integrate.simple(g, xs, xr.l, ppm.padding = ppm.padding)
    cbind(g, ints)
  })
  if (sum.dups > 0) warning(paste(sum.dups, "peaks were removed as duplicates. Merging not implemented."))
  
  pt = do.call("rbind", pt.l)
  xs@peaks = pt[order(pt[,"sample"]),]
  xs = buildGroups(xs, xs@peaks[,"new.gidx"])
  
  return(xs)
}


integrate.simple = function(g, xs, xr.l, ppm.padding = 1) {
  x = xs@peaks[g[,"pn"],,drop=F]
  
  mzrange.g = c(
    min(x[,"mzmin"], na.rm=T) - x[1,"mzmin"] * ppm.padding / 1E6, 
    max(x[,"mzmax"], na.rm=T) + x[1,"mzmax"] * ppm.padding / 1E6
  )
  
  laply(seq(nrow(g)), function(i) {
    p = x[i,]
    pg = g[i,]
    
    if (is.na(p["mzmin"]))  {
      mzrange = mzrange.g
    } else {
      mzrange = c(
        p["mzmin"] - p["mzmin"] * ppm.padding / 1E6, 
        p["mzmax"] + p["mzmin"] * ppm.padding / 1E6
        )
    }
    
    scanrange = as.numeric(floor(c(pg["sc.min"], pg["sc.max"])))
    scanrange[scanrange < 1] = 1
    maxscan = length(xr.l[[pg["sample"]]]@scantime)
    scanrange[scanrange > maxscan] = maxscan 
    eic = rawEIC(xr.l[[pg["sample"]]], mzrange = mzrange, scanrange = scanrange)
    
    # Recalc mz
    scans = llply(scanrange[1]:scanrange[2], function(x) getScan(xr.l[[pg["sample"]]], x, mzrange))
    scan.mat = do.call("rbind", scans)
    min = quantile(scan.mat[,"intensity"], .5)
    mz = mean(scan.mat[scan.mat[,"intensity"] > c(min), "mz"])
    
    rts = xs@rt$corrected[[pg[["sample"]]]][eic$scan]
    
    weightedi = diff(rts) * sapply(
      2:length(eic$intensity), function(x){
        mean(eic$intensity[c(x - 1, x)])
      })
    into = sum(weightedi)
    rtwmean = sum(weightedi * rts[-1]) / into
    
    as.matrix(data.frame(
      maxo = max(eic$intensity),
      into = into,
      mzmin.fill = mzrange[[1]],
      mzmax.fill = mzrange[[2]],
      rtmin = min(rts),
      rtmax = max(rts),
      rt = rtwmean,
      mz = mz
    ))
  })
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

plotGroup.xs = function(i, xs, xr.l) {
  g = xs@groupidx[[i]]
  ps = xs@peaks[g,,drop=F]
  
  scan.range = c(floor(min(ps[,"sc.min"], na.rm=T)), ceiling(max(ps[,"sc.max"], na.rm=T)))
  eic.mat = matrix(numeric(), ncol = nrow(ps), nrow = scan.range[2]-scan.range[1]+1)
  
  for (j in seq(ps[,1])) {
    mz.range = unname(c(ps[j,"mzmin.fill"], ps[j,"mzmax.fill"]))
    eic.mat[,j] = rawEIC(xr.l[[ps[j,"sample"]]], mzrange = mz.range, scanrange = scan.range)$intensity 
  }
  
  plot.me = melt(eic.mat)
  colnames(plot.me) = c("scan", "peak", "int")
  plot.me[,"peak"] = factor(plot.me[,"peak"])
  
  annos = data.frame(ps)
  annos[,"n"] = seq(annos[,1])
  colnames(annos)[1:4] = c("peak","scwm", "scmax", "scmin")
  annos[,2:4] = annos[,2:4] - min(scan.range)
  
  ggplot() + 
    geom_path(data = plot.me, aes(y=int, x= scan, colour=peak)) + 
    geom_point(data=annos, mapping=aes(y=0, x=scmax), colour="#000000") + 
    geom_point(data=annos, mapping=aes(y=0, x=scmin), colour="#000000") +
    geom_vline(data=annos, mapping=aes(xintercept=scwm), colour="#000000", alpha=0.3)+
    facet_grid(peak ~ ., scales="free_y")
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

add.raw.sc = function(xs) {
  rt.scans = matrix(NA, ncol=3, nrow=nrow(xs@peaks), dimnames=list(NULL, c("sc", "scmin", "scmax")))
  
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