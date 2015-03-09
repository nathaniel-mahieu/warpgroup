group.warpgroup = function(
  xs,
  xr.l=NULL,
  sc.drift.lim, 
  ppm.lim, 
  sc.lim, 
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
      .packages=c("warpgroup", "dtw", "igraph"),
      .export = c("xs", "xr.l")
      ), 
    xs, xr.l,
    
    .fun = function(i, xs, xr.l) {
      cat("\r", i)
      tryCatch(
        xcms.warpgroup(
          xs = xs, 
          g.i = i, 
          xr.l = xr.l,
          sc.drift.lim=sc.drift.lim, 
          ppm.lim = ppm.lim, 
          sc.lim = sc.lim
          ),
        error=paste,
        finally = return
        )
      }
    )
  
  if (output.groups) return(groups)
  return(warpgroupsToXs(xs, groups, xr.l)) # Merge back into traditional xcms workflow
}

warpgroupsToXs = function(xs, groups, xr.l) {
  cat("Converting warpgroups to xcmsSet.\nNote: The xcmsSet returned by this function does not need fillpeaks().\nCaution: diffreport() performs further processing on the peak groups before reporting statistics. Specifically it discards overlapping groups. This could remove groups which describe different portions of a peak found by the warpgrouping yet overlap.  If this behavior is not desired statstics can easily be performed on the raw warpgroup data retrieved by setting output.groups=T.\n")
  
  group.l = unlist(groups, F)
  bad.gs = which(!sapply(group.l, is.matrix))
  for (bg in rev(bad.gs)) group.l[[bg]] = NULL
  n.samps = nrow(group.l[[1]])
  
  pt = do.call("rbind", group.l)
  pt = cbind(pt, group = unlist(lapply(seq(group.l), rep, times=n.samps)))
  
  ints = llply(group.l, .progress="text", integrate.simple, xs, xr.l, ppm.padding=1)
  
  pt = cbind(pt, do.call("rbind", ints))
  pt = pt[order(pt[,"sample"]),]
  xs@peaks = pt
  
  xs = buildGroups(xs, xs@peaks[,"group"])
  
  return(xs)
}

xcms.warpgroup = function(
  xs, 
  g.i, 
  xr.l,
  sc.drift.lim, 
  ppm.lim, 
  sc.lim
  ) { #Handles xcms setup
  
  xs = add.raw.sc(xs)
  
  ps = xs@peaks[xs@groupidx[[g.i]],,drop=F] 
  ps = cbind(
    ps,
    n = 1:nrow(ps),
    sc.wmean = rep(NA, nrow(ps))
    )
 
  scans = c(
    min(ps[,"sc.min"])-sc.drift.lim, 
    max(ps[,"sc.max"])+sc.drift.lim
    )
  scans[scans <= 0] = 1
  start.scan = min(scans)-1
  
  #Per peak eic matrix
  eic.mat = matrix(integer(1), ncol=nrow(ps), nrow=length(scans[1]:scans[2]))
  for (i in seq_along(ps[,1])) {
    eic = rawEIC(xr.l[[ps[i,"sample"]]], mzrange = c(ps[i,"mzmin"], ps[i,"mzmax"]), scanrange = scans)$intensity 
    eic.mat[1:length(eic),i] = eic
  }
  
  #Find weighted mean of centWave bounds to help define peak groups
  sc.wmean = laply(seq_along(eic.mat[1,]), function(i) {
    eic = eic.mat[,i]
    peakrange = ps[i,"sc.min"]:ps[i,"sc.max"] - start.scan
    sum(eic[peakrange] * peakrange)/sum(eic[peakrange]) + start.scan
  })
  ps[,"sc.wmean"] = c(sc.wmean)
  
  
  #Per sample, summed EIC matrix for alignment
  mzmin = min(ps[,"mzmin"] - ps[,"mzmin"] * ppm.lim/1E6)
  mzmax = max(ps[,"mzmax"] + ps[,"mzmax"] * ppm.lim/1E6)
  
  eic.mat.s = matrix(integer(1), ncol=length(xr.l), nrow=length(scans[1]:scans[2]))
  for (i in seq_along(xr.l)) { 
    eic = rawEIC(xr.l[[i]], mzrange = c(mzmin, mzmax), scanrange = scans)$intensity
    eic.mat.s[1:length(eic),i] = eic
  }

  
  groups = warpgroup(
    ps = ps[,c("mz", "sc", "sc.wmean", "sc.min", "sc.max", "sample"),drop=F], # mz, sc, sc.wmean, scmin, scmax, sample
    eic.mat.s = eic.mat.s,
    start.scan = start.scan,
    
    sc.drift.lim,
    ppm.lim,
    sc.lim
  )
  
  groupidx = rep(g.i, ncol(eic.mat.s))
  pn = data.frame(n = seq(xs@groupidx[[g.i]]), pn = xs@groupidx[[g.i]])
  
  llply(groups, function(g) {
    pn = merge(g, pn, by="n", all.x = T)[,"pn"]
    cbind(g, groupidx, pn)
    })
}



integrate.simple = function(g, xs, xr.l, ppm.padding = 1) {
  x = xs@peaks[g[,"pn"],,drop=F]
  
  mzrange.g = c(
    mean(x[,"mzmin"], na.rm=T) - x[1,"mzmin"] * ppm.padding / 1E6, 
    mean(x[,"mzmax"], na.rm=T) + x[1,"mzmax"] * ppm.padding / 1E6
  )
  
  laply(seq(nrow(g)), function(i) {
    p = x[i,]
    pg = g[i,]
    
    if (is.na(p["mzmin"]))  {
      mzrange = mzrange.g
    } else {
      mzrange = c(p["mzmin"], p["mzmax"])
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


plotGroup = function(i, groups, integrates, xs, xr.l) {
  g = groups[[i]]
  int = integrates[[i]]
  ps = xs@peaks[na.omit(g[,"centWave.pn"]),,drop=F]
  
  scan.range = c(min(g[,"sc.min"], na.rm=T), max(g[,"sc.max"], na.rm=T))
  eic.mat = matrix(numeric(), ncol = nrow(g), nrow = scan.range[2]-scan.range[1])
  
  for (j in seq(g[,1])) {
    if (!is.na(g[j,"mzmin"])) {
      mz.range = c(g[j,"mzmin"], g[j,"mzmax"])
    } else {
      mz.range = c(int[j,"mzmin.fill"], int[j,"mzmax.fill"])
    }
    
    eic.mat[,j] = rawEIC(xr.l[[g[j,"samp"]]], mzrange = mz.range, scanrange = scan.range)$intensity 
  }
  
  plot.me = melt(eic.mat)
  colnames(plot.me) = c("scan", "peak", "int")
  plot.me[,"peak"] = factor(plot.me[,"peak"])
  
  annos = dcast(melt(int))
  
  ggplot(plot.me, aes(y=int, x= scan, colour=peak)) + 
    geom_path() + 
    facet_grid(peak ~ ., scales="free_y")
  
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