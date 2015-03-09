group.stats = function(xs, groups, integrates, samples.a, samples.b) {
  laply(integrates, function(x) {
    as = subset(x, samp %in% samples.a)
    bs = subset(x, samp %in% samples.b)
    
    data.frame(
      mean.a = mean(as[,"into"]),
      sd.a = sd(as[,"into"]),
      mean.b = mean(bs[,"into"]),
      sd.b = sd(bs[,"into"]),
      pval = t.test(as[,"into"], bs[,"into"])$p.value,
      fold = mean(as[,"into"])/mean(bs[,"into"]),
      n.centWave.a = sum(!is.na(as[,"centWave.pn"])),
      n.centWave.b = sum(!is.na(bs[,"centWave.pn"]))
      )
    })
}

format.output = function(groups, stats, integrates) {
  
  laply(seq(groups), function(i) {
    data.frame(
      mz = mean(groups[[i]][,"mz"], na.rm=T),
      group = i,
      integrates[i,],
      sc.wmean = mean(stats[[i]][,"sc.wmean"], na.rm=T)
      )
    })
  
}