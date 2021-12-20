# Popsom
### An Efficient Implementation of Kohonen's Self-Organizing Maps (SOMs) with Starburst Visualizations

![](https://raw.githubusercontent.com/lutzhamel/popsom/master/map.png)

Kohonen's self-organizing maps with a number of distinguishing features:

1. [An efficient, single threaded, stochastic training algorithm](https://doi.org/10.1007/978-3-030-01057-7_60) inspired by ideas from tensor algebra.  No special accelerator hardware required.

2. Automatic centroid detection and [visualization using starbursts](http://citeseerx.ist.psu.edu/viewdoc/citations;jsessionid=ADD8A6B07C6E90AC162CC1BCA57E996E?doi=10.1.1.217.6666).

3. Two models of the data: (a) a self organizing map model, (b) a centroid based clustering model.

4. A number of easily accessible [quality metrics for the self organizing map](https://doi.org/10.1007/978-3-319-28518-4_4) and the centroid based cluster model.

For an example usage case check [here](https://www.kaggle.com/lutzhamel/customer-segmentation-with-soms).

For the Popsom description as part of CRAN check [here](https://CRAN.R-project.org/package=popsom).
