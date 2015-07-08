# COntinued from X:\Nathaniel Mahieu\Q Exactive Plus\2NM109A Anno\3\creDBle_processing_1435188608\warpgroup

xs.bak = xs



#10140 - bounds messed up
#14720 - wrong peaks
#29591 - noise

xs = xs.bak
xs@groupidx = xs@groupidx[i]

xs = add.raw.sc(xs)

myiter = iter.gwparams(xs, xr.l, sc.max.drift = 70, ppm.max.drift = 3)

params = nextElem(myiter)


groups = warpgroup(
  ps = params$ps, # sc, scmin, scmax, sample
  eic.mat.s = params$eic.mat.s,

  sc.aligned.lim = 10,
  sc.aligned.factor = 1,
  detailed.groupinfo=T,
  min.peaks = 1,
  tw="dtw"
)
groups = llply(groups, function(x) {
  x[,c("sc", "scmin", "scmax")] = x[,c("sc", "scmin", "scmax")] + params$start.scan
  
  cbind(x, pn = params$gidx[x[,"n"]])
})
groups

  wg = groups
  wg = lapply(wg, function(x) { x[!duplicated(x[,"sample"]),] })
  plots = lapply(wg, plot.warpgroup, xs2, xr.l, sc.pad = 100, mz.pad = 0.001)
  do.call(grid.arrange, plots)


ps = params$ps
eic.mat.s = params$eic.mat.s
sc.aligned.lim = 10
sc.aligned.factor = 1
detailed.groupinfo=T
min.peaks = 1
tw="dtw"
pct.pad = 0.1

