# Common Questions
- **Where is warp.consistency?** This parameter is computed when the flag detailed.groupinfo is set to true. Eg. `group.warpgroup(detailed.groupinfo=T)`

- **Warpgroup never finishes running.** Always test warpgroup on your specific settings and data with a subset to make sure it will complete in a reasonable amount of time.  Due to poor scaling it is possible your analysis will never complete. Ex.
```r
xs2 = xs
xs2@groupidx = sample(xs2@groupidx, 10)
start = Sys.time()
group.warpgroup(xs2, ...)
end = Sys.time()
paste("Full dataset will take", length(xs@groupidx)/length(xs2@groupidx), "times as long to process.")
paste("Subset took", end - start)
```

- **XCMS Function xxxxx returns an error** Warpgroup returns the output of the algorithm and packs it in an XCMS data structure.  Warpgroup does not ensure that everything XCMS expects is true of the data.  For example, if no ion intensity was detected for one peakgroup in a sample Warpgroup does not assign a peak m/z or retention time to that sample (there was none). Most known errors are due to NA values in the peak table.  Try removing all peaks with NA values and rebuilding the groupvals with something like below.
```r
missing.values = which(is.na(xs@peaks[,c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into")]), arr.ind = T)
for (i in seq(nrow(missing.values))) {
  ps = xs@groupidx[[xs@peaks[i,"new.gidx"]]]
  xs@peaks[i, missing.values[i,]] = mean(xs@peaks[ps, missing.values[i,]], na.rm=T)
}

xs = warpgroup::refreshGroups(xs)
```

# Parallel Computing
Warpgroup supports parallel computation of a job via the foreach backend.  Simply mount any foreach compatible backend and warpgroup will automatically use that to distribute the job.  doParallel is great for local cores and doRedis is great for distributing the job among multiple machines. For parallel computing something like the following example code should be run before starting group.warpgroup. 

## doParallel
```r
library(doParallel)
cl = makeCluster(detectCores() - 1)
registerDoParallel(cl)
```

## doRedis
Master node
```r
library(doRedis)
registerDoRedis("worklist.name", "IP.OF.REDIS.SERVER")
```

Processing node(s)
```r
library(doRedis)
startLocalWorkers(n, "worklist.name", "IP.OF.REDIS.SERVER")
```

# Thank You
Thank you to the following for beta testing and bug report assistance.
- Kevin Cho
- Robert Sander Jansen