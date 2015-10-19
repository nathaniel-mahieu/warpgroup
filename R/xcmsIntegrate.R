#' A wrapper for integration of warpgroups and packaging as an xcmsSet object.
#' 
#' \code{warpgroupsToXs} returns an xcmsSet with the original groups and peaks replaced by those described in \code{groups}.
#' 
#' This function simply iterates over each group, integrates the appropriate region from xr.l and packages the results as an xcmsSet.

#' @param xs An xcmsSet object.
#' @param xr.l A list contanining xcmsRaw objects for each sample in the xcmsSet object in order.
#' @param groups The list of warpgroups as output by \code{group.warpgroup(output.groups=TRUE)}
#' @param ppm.padding Float. An expansion of the mass window for integration.
#' @param min.ppm.width Float. Any mass windows narrower than this are expanded to this width.
#' 
#' @return An xcmsSet with groups and peaks replaced by warpgroup generated ones.
#' 
#' @seealso See \url{https://github.com/nathaniel-mahieu/warpgroup} for examples.

warpgroupsToXs = function(xs, groups, xr.l, ppm.padding=0.1, min.ppm.width = 0) {
  cat("Converting warpgroups to xcmsSet.\nNote: The xcmsSet returned by this function does not need fillpeaks().\n")
    
  group.l = unlist(groups, F)
  bad.gs = unique(c(
    which(!sapply(group.l, is.matrix) | is.null(sapply(group.l, nrow)))
  ))
  if (length(bad.gs) > 0) warning(paste(length(bad.gs), "groups were removed due to errors."))
  for (bg in rev(bad.gs)) group.l[[bg]] = NULL
  
  pt.l = foreach(
    params = iter.integrateparams(group.l, xs, xr.l, ppm.padding, min.ppm.width),
    .errorhandling = "pass",
    .noexport = c("xr.l"),
    .inorder=T
  ) %do% {
    ints = integrate.simple(params)
    
    cbind(params[[2]], ints)
  }
  
  pt.l = lapply(seq(pt.l), function(i) cbind(pt.l[[i]], new.gidx = i))
  pt = do.call("rbind", pt.l)
  
  rts = matrix(numeric(), ncol = 2, nrow = nrow(pt), dimnames = list(NULL, c("rtmin", "rtmax")))
  rt = rep(numeric(), nrow(pt))
  xs@peaks = pt[order(pt[,"sample"]),]
  xs = buildGroups(xs, xs@peaks[,"new.gidx"])
  
  return(xs)
}

iter.integrateparams = function(group.l, xs, xr.l, ppm.padding, min.ppm.width = 0) {
  it <- iter(group.l)
  
  nextEl = function() {
    g <- nextElem(it) 
    
    g = g[!duplicated(g[,"sample"]),,drop=F]
    
    x = xs@peaks[g[,"pn"],,drop=F]
    
    mzrange.g = c(
      min(x[,"mzmin"], na.rm=T) - max(x[,"mzmin"], na.rm=T) * ppm.padding / 1E6, 
      max(x[,"mzmax"], na.rm=T) + max(x[,"mzmax"], na.rm=T) * ppm.padding / 1E6
    )
    
    mmzr = mean(mzrange.g)
    if (diff(mzrange.g)/mmzr * 1E6 < min.ppm.width) {
      mzrange.g = c(
        mmzr - mmzr * min.ppm.width/2/1E6,
        mmzr + mmzr * min.ppm.width/2/1E6
        )
    }
    
    pdata = foreach(i=seq(nrow(g))) %do% {
      p = x[i,]
      pg = g[i,]
      xr = xr.l[[pg["sample"]]]
      
      mzrange = mzrange.g
      if (!is.na(p["mzmin"])) mzrange = unname(c(p["mzmin"] - p["mzmin"] * ppm.padding / 1E6, p["mzmax"] + p["mzmin"] * ppm.padding / 1E6))
      ppm.width = abs(diff(mzrange)/mean(mzrange) * 1E6)
      if (ppm.width < min.ppm.width) mzrange = c(mzrange[1] - (abs(ppm.width - min.ppm.width)/2 * mzrange[1] / 1E6), mzrange[2] + (abs(ppm.width - min.ppm.width)/2 * mzrange[2] / 1E6))
      
      scanrange = as.numeric(c(
        which.min(abs(xr@scantime - floor(pg["rtmin.raw"]))), 
        which.min(abs(xr@scantime - ceiling(pg["rtmax.raw"])))
      ))
      scanrange[scanrange < 1] = 1
      maxscan = length(xr@scantime)
      scanrange[scanrange > maxscan] = maxscan 
      eic = rawEIC(xr, mzrange = mzrange, scanrange = scanrange)      
      eic = cbind(do.call("cbind", eic), rt=xs@rt$corrected[[pg["sample"]]][eic$scan])
      
      scans = foreach(i=scanrange[1]:scanrange[2]) %do% getScan(xr, i, mzrange)
      scanmat = do.call("rbind", scans)
      
      list(
        eic = eic,
        scanmat = scanmat,
        mzrange = mzrange,
        sc = scanrange
      )
    }
    
    list(
      pdata,
      g
      )
    
  }
  
  obj <- list(nextElem=nextEl)
  class(obj) <- c('ixts', 'abstractiter', 'iter')
  obj
}

integrate.simple = function(params) {
  int.mat = matrix(numeric(), nrow=length(params[[1]]), ncol=12, dimnames=list(NULL, c("mz", "rt", "rt.half", "rtmin", "rtmax", "into", "maxo", "mzmin", "mzmax", "mino", "scmin.raw", "scmax.raw")))
  
  for (i in seq(params[[1]])) {
    scan.mat = params[[1]][[i]]$scanmat
    eic = params[[1]][[i]]$eic
    
    min = quantile(scan.mat[,"intensity"], .5) * .9
    mz = mean(scan.mat[scan.mat[,"intensity"] > c(min), "mz"])
    
    weightedi = diff(eic[,"rt"]) * sapply(
      2:length(eic[,"intensity"]), function(x){
        mean(eic[,"intensity"][c(x - 1, x)])
      })
    into = sum(weightedi)
    rtwmean = sum(weightedi * eic[,"rt"][-1]) / into
    
    rt.half = mean(c(eic[,"rt"][1], eic[,"rt"][length(eic[,"rt"])]))
    
    int.mat[i,] = c(
      mz,
      rtwmean,
      rt.half,
      min(eic[,"rt"]),
      max(eic[,"rt"]),
      into,
      max(eic[,"intensity"]),
      params[[1]][[i]]$mzrange[[1]],
      params[[1]][[i]]$mzrange[[2]],
      min(eic[,"intensity"]),
      params[[1]][[i]]$sc[1],
      params[[1]][[i]]$sc[2]
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


refreshGroups = function(xs, keepall=F) {
  gs = sapply(seq(xs@groupidx), function(i) { rep(i, length(xs@groupidx[[i]])) })
  xs2 = xs; xs2@peaks = xs@peaks[unlist(xs@groupidx),]
  
  xs3 = buildGroups(xs2, gs)
  
  if (keepall) xs3@peaks = rbind(xs3@peaks, xs@peaks[-unlist(xs@groupidx),])
  
  xs3
}