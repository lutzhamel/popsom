# Popsom
### An Efficient Implementation of Kohonen's Self-Organizing Maps (SOMs) with Starburst Visualizations

![](https://raw.githubusercontent.com/lutzhamel/popsom/master/map.png)

Kohonen's self-organizing maps with a number of distinguishing features:

1. An efficient, single threaded, stochastic training algorithm inspired by ideas from tensor algebra.  Provides significant speedups over traditional single-threaded training algorithms. No special accelerator hardware required (see <doi:10.1007/978-3-030-01057-7_60>).

2. Automatic centroid detection and visualization using starbursts.

3. Two models of the data: (a) a self organizing map model, (b) a centroid based clustering model.

4. A number of easily accessible quality metrics for the self organizing map and the centroid based cluster model (see <doi:10.1007/978-3-319-28518-4_4>).

For an example usage case check [here](https://www.kaggle.com/lutzhamel/customer-segmentation-with-soms).

For the Popsom description as part of CRAN check [here](https://CRAN.R-project.org/package=popsom).
