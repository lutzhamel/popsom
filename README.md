# Popsom
### A Very Efficient Implementation of Kohonen's Self-Organizing Maps (SOMs) with Starburst Visualizations

![](https://raw.githubusercontent.com/lutzhamel/popsom/master/map.png)

POPSOM is a package implementing Kohonen's self-organizing maps with a number of
distinguishing features:
1. A very efficient, single threaded, stochastic training algorithm
based on ideas from tensor algebra.  Up to 60x faster than traditional single-threaded
training algorithms. No special accelerator hardware required.

2. Automatic centroid detection and visualization using starbursts.

3.  Maintains two models of the data: (a) a self-organizing map model, (b) a centroid based clustering model.

4.  Provides a number of easily accessible quality metrics for the self-organizing map and
the centroid based cluster model.

For an example usage case check [here](https://www.kaggle.com/lutzhamel/customer-segmentation-with-soms).

For the Popsom description as part of CRAN check [here](https://CRAN.R-project.org/package=popsom).
