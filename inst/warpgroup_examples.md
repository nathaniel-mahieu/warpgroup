# Warpgroup Examples

This is a set of example Warpgroup uses, each illustrating a different feature of the algorithm.




```r
library(warpgroup)
library(ggplot2)
```


#Data Format

```r
data(example_1_3)

head(eic.mat)
```

```
##          [,1]     [,2]      [,3]     [,4]     [,5] [,6] [,7] [,8]
## [1,]     0.00     0.00  9305.719 21281.18 11959.37    0    0    0
## [2,] 10862.48 11634.03     0.000 21920.07     0.00    0    0    0
## [3,]     0.00     0.00 12736.042 14096.86     0.00    0    0    0
## [4,]     0.00     0.00     0.000 34584.13     0.00    0    0    0
## [5,]     0.00     0.00     0.000 19327.64     0.00    0    0    0
## [6,] 21487.79     0.00     0.000 40626.04     0.00    0    0    0
```
Extracted Ion Chromatograms representative of the supplied peak bounds are supplied as a matrix, rows representing scans and columns representing samples.


```r
peak.bounds
```

```
##      sc scmin scmax sample
## [1,] 32    15    48      1
## [2,] 33    16    48      2
## [3,] 30    15    46      3
## [4,] 34    11    50      4
## [5,] 37    21    55      5
## [6,] 36    19    53      6
## [7,] 37    20    55      7
## [8,] 33    16    50      8
```
Peak bounds detected in a prerequisite step are supplied as a matrix.  Columns scmin, scmax, and sc are integers refering to rows of the EIC matrix. Column sample refers to the column in the EIC matrix.

#1. Consensus Peak Bound Determination

```r
data(example_1_3)
plot_peaks_bounds(eic.mat, peak.bounds)
```

<img src="figure/bound_variation_start-1.png" title="plot of chunk bound_variation_start" alt="plot of chunk bound_variation_start" style="display: block; margin: auto;" />

Our initial EICs and bounds all agree well (they were previously warpgrouped).  Lets add some variance to the peak bounds to simulate real world peak detection.


```r
bad.bounds = aperm(apply(peak.bounds, 1, function(x) {
  x[2] = x[2] + floor(runif(1, 5,11)) * sample(c(-1,1),1)
  x[3] = x[3] + floor(runif(1, 5,11)) * sample(c(-1,1),1)
  x
  }))
plot_peaks_bounds(eic.mat, bad.bounds)
```

<img src="figure/bound_variation_challenge-1.png" title="plot of chunk bound_variation_challenge" alt="plot of chunk bound_variation_challenge" style="display: block; margin: auto;" />

We see the peak bounds are not self-consistent, each set of bounds describing a different peak region.  This artifically increases our measurement variance and decreases our statstical power.

Lets see if we can fix this.


```r
wg.bounds = warpgroup(bad.bounds, eic.mat, sc.aligned.lim = 6)[[1]]
plot_peaks_bounds(eic.mat, wg.bounds)
```

<img src="figure/bound_variation_fixed-1.png" title="plot of chunk bound_variation_fixed" alt="plot of chunk bound_variation_fixed" style="display: block; margin: auto;" />

After submitting the EIC traces and detected peak bounds to warpgrouping we see that all the samples have similar integration regions detected.

#2. Peak Subregion Detection
Still need to find a good example of this

#3. Determination of Undetected Peak Bounds (Fillpeaks)

```r
data(example_1_3)
plot_peaks_bounds(eic.mat, peak.bounds)
```

<img src="figure/missing_bounds_start-1.png" title="plot of chunk missing_bounds_start" alt="plot of chunk missing_bounds_start" style="display: block; margin: auto;" />

Again our inital EICs all look good.  Lets delete all the peaks but one and see how warpgroup handles this.


```r
bad.bounds = peak.bounds[1,,drop=F]
plot_peaks_bounds(eic.mat, bad.bounds)
```

<img src="figure/missing_bounds_deleted-1.png" title="plot of chunk missing_bounds_deleted" alt="plot of chunk missing_bounds_deleted" style="display: block; margin: auto;" />

We see only one sample had a peak detected in this region. Traditionally a large, indiscriminate region of the chromatogram is integrated to fill these missing values.  The warpgroup solution follows.


```r
wg.bounds = warpgroup(bad.bounds, eic.mat, sc.aligned.lim = 6)[[1]]
plot_peaks_bounds(eic.mat, wg.bounds)
```

<img src="figure/missing_bounds_fixed-1.png" title="plot of chunk missing_bounds_fixed" alt="plot of chunk missing_bounds_fixed" style="display: block; margin: auto;" />

After submitting the EIC traces and detected peak bound to warpgrouping we see that all the samples now have similar integration regions.


#4. Grouping Peaks Which Deviate From Global Retention Time Correction

```r
data(example_4)

plot_peaks_bounds(eic.mat, peak.bounds)
```

<img src="figure/mixed_group_start-1.png" title="plot of chunk mixed_group_start" alt="plot of chunk mixed_group_start" style="display: block; margin: auto;" />

This time there are two peaks per sample in this rough group.  These have been warpgrouped, but in real life samples each peak's drift will deviate somewhat from the global retention time shift.  Lets add some extreme drift.


```r
eic.mat.bad = array(numeric(), dim=dim(eic.mat))
peak.bounds.bad = peak.bounds

for (r in seq(ncol(eic.mat))) {
  shift = floor(runif(1,0,470))
  eic.mat.bad[,r] = c(
    eic.mat[shift:nrow(eic.mat),r], 
    rep(0,times=shift-1)
    )
  
  peak.bounds.bad[peak.bounds.bad[,"sample"] == r,1:3] = peak.bounds[peak.bounds[,"sample"] == r,1:3] - shift
}

plot_peaks_bounds(eic.mat.bad, peak.bounds.bad)
```

<img src="figure/mixed_group_drifting-1.png" title="plot of chunk mixed_group_drifting" alt="plot of chunk mixed_group_drifting" style="display: block; margin: auto;" />

Here is a case that current algorithms would have no chance at grouping. The inter-sample drift is greater than the separation between the peaks to be grouped.  Lets see how warpgrouping does.


```r
wg.bounds = warpgroup(peak.bounds.bad, eic.mat.bad, sc.aligned.lim = 6)

for (g in wg.bounds) print(plot_peaks_bounds(eic.mat.bad, g))
```

<img src="figure/mixed_group_fixed-1.png" title="plot of chunk mixed_group_fixed" alt="plot of chunk mixed_group_fixed" style="display: block; margin: auto;" /><img src="figure/mixed_group_fixed-2.png" title="plot of chunk mixed_group_fixed" alt="plot of chunk mixed_group_fixed" style="display: block; margin: auto;" />

We see that warpgroup has correctly grouped each peak into the proper, distinct group.

#5. An Extreme Example
This is an extreme example, data this unreliable probably shouldn't be trusted... but it provides a nice challenge.


```r
data(example_5)

plot_peaks_bounds(eic.mat, peak.bounds)
```

<img src="figure/extreme_example-1.png" title="plot of chunk extreme_example" alt="plot of chunk extreme_example" style="display: block; margin: auto;" />

We can clearly see two peaks in most samples.  There is a large retention time drift.  There is also a varying degree of merging between the two peaks.  In some samples two distinct peaks were detected, in others a single peak was detected.  Lets see how warpgroup handles this.


```r
wg.bounds = warpgroup(peak.bounds, eic.mat, sc.aligned.lim = 8)

for (g in wg.bounds) print(plot_peaks_bounds(eic.mat, g))
```

<img src="figure/extreme_example_fixed-1.png" title="plot of chunk extreme_example_fixed" alt="plot of chunk extreme_example_fixed" style="display: block; margin: auto;" /><img src="figure/extreme_example_fixed-2.png" title="plot of chunk extreme_example_fixed" alt="plot of chunk extreme_example_fixed" style="display: block; margin: auto;" /><img src="figure/extreme_example_fixed-3.png" title="plot of chunk extreme_example_fixed" alt="plot of chunk extreme_example_fixed" style="display: block; margin: auto;" />

Warpgroup grouped the peaks into three distinct groups, each describing a different chromatographic region.

