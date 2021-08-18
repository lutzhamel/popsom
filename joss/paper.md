---
title: 'popsom: An Efficient Implementation of Self-Organizing Maps with Starburst Visualizations for R'
tags:
  - R
  - machine learning
  - unsupervised learning
  - neural networks
authors:
  - name: Lutz Hamel
    orcid: 0000-0002-2302-3358
    affiliation: 1
affiliations:
 - name: Department of Computer Science and Statistics, University of Rhode Island, Kingston, RI 02881
   index: 1
date: 7 July 2021
bibliography: paper.bib
---

# Background

The self-organizing map (SOM) was developed by Teuvo Kohonen for data exploration and
data visualization in the 1980s [@kohonen2001self]. Specifically, the self-organizing map
is an unsupervised artificial neural network designed for identifying and visualizing clusters in multi-dimensional
spaces.  What makes the self-organizing map so attractive are the intuitive mathematical
underpinnings and the straightforward visualization of computational
results [@ultsch1990self]. Self-organizing maps have been applied in virtually every
scientific discipline where some sort of data exploration or analysis is
necessary, e.g., [@liakos2018machine; @mathys2019single; @matic2018interpreting;
@miller1996star].  A number of R-packages exist that implement self-organizing maps, e.g.,
[@kohonen2018jss; @som2016R].


# Statement of need

Training a self-organizing map is time consuming. Here
we introduce an R-package called `popsom` [@popsom2021R] that implements a training algorithm
for self-organizing maps inspired by
tensor algebra [@hamel2018vsom] which provides significant speedups over
traditional implementations of the training algorithm such as the `som` package by Yan
[@som2016R].  These speedups enable researchers to look at much larger datasets or to improve
their throughput with a dataset of a given size.

Visualization of a trained map is at the core of using self-organizing maps.  The
`popsom` package improves on the standard u-matrix visualization
for self-organizing maps [@ultsch1990self] by superimposing starbursts in order to
highlight cluster structures [@hamel2011improved].

# Description

At a slightly more detailed level, our `popsom` package implements
Kohonen's self-organizing maps with a number of distinguishing features:

1. A very efficient, single threaded, stochastic training algorithm based on ideas from tensor algebra providing significant speedups over traditional single-threaded implementations
without the need for special accelerator hardware. The speedups result from the fact that the
vector and matrix structures exposed by our algorithm map neatly into
vector and matrix operations available on today's CPUs.  Our Fortran 90 implementation
insures that these vector and matrix operations are mapped onto the hardware as efficiently as
possible [@hamel2018vsom].

2. Automatic centroid visualization and detection using starbursts [@hamel2011improved]. Not
only does `popsom` display clusters and centroids on the map using starbursts but it
also computes a cluster model based on the starbursts.  

3. The `popsom` package maintains two models of the given training data: (a) a
self-organizing map model where elements of the map model are available to
the user for analysis, and (b) a centroid based clustering model similar to a k-means
model where centroid and cluster information is available to the user.  Having these
two perspectives of a dataset is often helpful during a data analysis.

4. The package provides a number of easily accessible quality metrics for the self-organizing map and the centroid based cluster models [@hamel2016som; @tatoian2018self]. In particular, the package computes the `convergence` of a map which is a linear combination of the variance captured and the topographic fidelity of the map. A value close to one of this metric indicates a converged map. Furthermore, `popsom` also computes the `separation` of the clusters
in a model.  In general, a value close to one here indicates well separated clusters.

# Usage

`popsom` is available on [CRAN](https://CRAN.R-project.org/package=popsom) and can be installed and loaded into an R session using,
```
> install.packages("popsom")
> library(popsom)
```
Binary packages for `popsom` are available from CRAN for macOS, Linux, and Windows.  If you are on a system
that is not supported by CRAN, you can download and compile the package from
[GitHub](https://github.com/lutzhamel/popsom).

The following is a simple use-case for `popsom` exercising some of the functionality
it has to offer. We start by constructing a model,
```
> ## load a dataset
> data(iris)
>
> ## set data frame and labels
> df <- subset(iris,select=-Species)
> labels <- subset(iris,select=Species)
>
> ## build a self-organizing map
> m <- map.build(df,labels,xdim=15,ydim=10,train=100000,seed=42)
>
> ## compute a summary and display it
> map.summary(m)

Training Parameters:
  xdim ydim alpha train normalize seed instances columns
    15   10   0.3 1e+05     FALSE   42       150       4

Quality Assessments:
  convergence separation clusters
         0.94       0.93        4

>
```
The `map.summary` function gives us a quick snapshot of relevant model information. Perhaps
the most interesting thing here are the `Quality Assessments`. The `convergence` value is a linear
combination of the estimated topographic accuracy and the embedding accuracy.  The latter roughly
corresponds to the amount of training data variance the map models [@hamel2016som].  Convergence values of 0.9 and
better are considered hallmarks of high quality maps.  The `separation` value is computed with the
formula,
```
1 - wcss/bcss
```
where `wcss` is the average within cluster sum of squares and `bcss` is the average between cluster sum
of squares.  This computation is a quick way of assessing
the quality of the computed clusters.  Here, a value close
to one usually indicates a good cluster model.

We can look at details of the centroids of the model,
```
> ## look at the centroids of the model
> lc <- m$unique.centroids
> ll <- m$centroid.labels
> for (i in 1:length(lc))
+ cat("(",lc[[i]]$x,",",lc[[i]]$y,") -> ",ll[[lc[[i]]$x,lc[[i]]$y]],"\n")
( 4 , 1 ) ->  versicolor
( 1 , 4 ) ->  virginica
( 9 , 1 ) ->  versicolor
( 15 , 10 ) ->  setosa
>
```
The map model maintains information about the centroids which we can access.  Here
we access that information in order to print out the coordinates of the centroids on the map
together with their assigned labels.  We can easily verify this information using
the starburst plots the package provides,
```
> ## display a starburst plot of the map model
> map.starburst(m)
```
The starburst visualization shown
in \autoref{fig:map}.  The centroids we extracted from the model earlier are easily
identified on this visualization.  The starbursts indicate the extent of clusters around
each centroid and the colors of the heatmap indicate the "tightness" of the clusters.
Hot, yellow usually indicates tight clusters whereas reddish colors indicate looser clusters or
borders of clusters.  Here we can see that the "versicolor" and "virginica" clusters are more similar
to each other and that the "setosa" cluster is separated from the other clusters by some
distance. This is easily verified with a scatter plot matrix of the iris dataset using
```
> colors <- c("#00AFBB", "#E7B800", "#FC4E07")  
> pairs(iris[,1:4], pch = 19,  cex = 0.5,
+       col = colors[iris$Species])
>
```
and shown in \autoref{fig:scatter}.  The blue cluster is the "setosa" cluster and it is easy to see that
it is well separated from the other clusters.

A more involved usage example can be found on
[Kaggle](https://www.kaggle.com/lutzhamel/self-organizing-maps-in-customer-segmentation).

![Starburst visualization of a self-organizing map.\label{fig:map}](map.png)

![Scatterplot matrix of the iris dataset.\label{fig:scatter}](scatter.png)

# Performance

In order to highlight the kind of performance gains our package typically provides we
included a simple
[benchmark script](https://github.com/lutzhamel/popsom/blob/master/package/performance/popsom-perf.r)
in our repository.
The benchmark compares the training algorithm performances of three packages: (1)  `popsom` [@popsom2021R],
(2) `som` [@som2016R], and (3) `kohonen` [@kohonen2018jss].
The script uses three real world datasets and one synthetic dataset.
We only consider training times of converged maps, that is, we only look at the training times of
maps that have high topographic accuracy
and capture most of the variance of the input data. In the terminology of our `popsom` package
that means we only consider training times of maps with a convergence value of 0.9 or better.
For the output below,
the benchmarks were run on an Intel Linux machine and the output was slightly edited for readability.
The benchmarks are roughly ordered in increasing complexity of the data.  We start with the `iris` dataset.

## The Iris Dataset

The `iris` dataset is available in R. A 15x10 map size and 100000 training iterations are sufficient to consistently produce
converged maps which the following summary confirms,
```
Training Parameters:
  xdim ydim alpha train normalize seed instances columns
    15   10   0.3 1e+05     FALSE   42       150       4

Quality Assessments:
  convergence separation clusters
         0.94       0.93        4

```
This is the simplest dataset in our benchmark with 150 instances and 4 columns.
Running the benchmarks gives us the following performance data where the time is reported in milliseconds,
```
  package        mean time          speedup
1  popsom           214.96                1
2     som          4092.20               19
3 kohonen         23535.20              109
```
Here we can see that `popsom` is 19 times faster than `som` and 109 times faster than `kohonen`.

## The Epil Dataset

The `epil` dataset is available in R through the `MASS` package.  Here again, a map size of 15x10 and
100000 training iterations are sufficient to produce converged maps.  The following is a summary
of building a converged map,
```
Training Parameters:
  xdim ydim alpha train normalize seed instances columns
    15   10   0.3 1e+05     FALSE   42       236       8

Quality Assessments:
  convergence separation clusters
         0.97       0.83        4
```
With 236 instances and 8 columns the dataset is slightly more complex than the `iris` dataset.
The results from the benchmark are,
```
  package        mean time         speedup
1  popsom           397.73               1
2     som          4740.05              12
3 kohonen         23466.17              59
```
Times are reported in milliseconds.  Here we can see that `popsom` is about 12 times faster than the `som`
package and about 60 times faster than the `kohonen` package.

## The Wines Dataset

The `wines` dataset is available through the `kohonen` package. Like in the previous two benchmarks,
a map size of 15x10 and
100000 training iterations are sufficient to consistently produce converged maps.  The following is a summary
of building a converged map,
```
Training Parameters:
  xdim ydim alpha train normalize seed instances columns
    15   10   0.3 1e+05     FALSE   42       177      13

Quality Assessments:
  convergence separation clusters
         0.96       0.93        4
```
What distinguishes this dataset from the others is the relative high dimensionality compared to the other
datasets.  The results for the benchmark are,
```
  package        mean time          speedup
1  popsom           631.33                1
2     som          5453.73                9
3 kohonen         23433.91               37
```
Times are reported in milliseconds. Here we measure a speedup of 9 compared to the `som` package and a speedup of close to 40 compared
to the `kohonen` package.

## The Synthetic Dataset

Our script constructs a synthetic dataset containing three clusters embedded in an 10-dimensional space. With
300 observations it represents our largest and most complex dataset. In order to obtain
converged maps consistently we had to pick a map size of 25x20 together with 1000000 training iterations.
The following is a summary of building a map with these parameters,
```
Training Parameters:
  xdim ydim alpha train normalize seed instances columns
    25   20   0.3 1e+06     FALSE   42       300      10

Quality Assessments:
  convergence separation clusters
            1          1        9
```
As one would expect from a synthetic dataset without any noise the convergence value
for this map is one: no topographic errors and full modeling of the underlying variance.
Even though the summary reports 9 clusters have been found it is easily verified visually through
the use of the starburst function that those 9 clusters are actually subclusters of three major clusters
as required by the design of the synthetic data.  The performance summary is,
```
  package        mean time          speedup
1  popsom            15.48                1
2     som           164.40               10
3 kohonen          2641.17              170
```
Here we can see that `popsom` is about an order of magnitude faster than the `som` package and a little bit
more than two orders of magnitude faster than the `kohonen` package.
Times are reported in seconds.

## Observations

With our package one can expect roughly an order of magnitude speedup over
the `som` package and speedups between one and two orders of magnitude over
the `kohonen` package for most datasets.  Given that these packages represent typical
implementations of single-threaded, stochastic training algorithms it is highly likely that
our speedup numbers also apply to other SOM implementations that implement
single-threaded, stochastic training of self-organizing maps.

As we saw in the benchmarks above, the exact speedup provided by our package
depends on the precise constellation of the training data as well as the map size
and training iterations to provide consistently converged maps.
Therefore it is difficult to come up with a precise characterization of our speedup
numbers.  However, it is evident that our algorithm is sensitive to the dimensionality
of the training data.  But even this characterization is not as straightforward as it
might seem since higher dimensional data often requires larger maps and more
training iterations in order to converge giving our algorithm additional opportunities
for optimization.

A more detailed yet older performance analysis of our `popsom` package appears in [@hamel2018vsom].

# Acknowledgements

We would like to thank the University of Rhode Island for financial support of
this project through their 'Project Completion Grants.' A special thanks to the
following people for their contributions to the project
in no particular order: Benjamin Ott, Gregory Breard,  Robert Tatoian,
Michael Eiger, and Vishakh Gopu.  We also would like to thank sthda.com for their
[wiki on scatterplots](http://www.sthda.com/english/wiki/scatter-plot-matrices-r-base-graphs).

# References
