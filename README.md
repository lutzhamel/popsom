# Popsom (Version 5.0)
![](https://raw.githubusercontent.com/lutzhamel/popsom/master/map.png)

Popsom is a R package for the creation and evaluation of self-organizing maps (SOMs).  Popsom version 5.0 supports the following features:
- Support for two models:
  1. A self-organizing map model
  2. A centroid based clustering model
- Quality measures available for both models
- Streamlined S3 based API
  - Easy access to the most important map and centroid data structures
- Powerful map visualization with centroid identification
- Extremely fast training algorithm based on ideas from tensor algebra

For an example usage case check [here](https://www.kaggle.com/lutzhamel/customer-segmentation-with-soms).

For the Popsom description as part of CRAN check [here](https://CRAN.R-project.org/package=popsom).

# Popsom (Version 5.1)

- Something got rattled with the S3 interface in R 4.x.  It now longer works the way it did in release 3.x.  Therefore, I took the S3 interface out because I want the package to work with both 3.x and 4.x installations.  Furthermore, the advantages of the S3 interface are incremental at best and I don't feel like debugging R internals.
- Implemented a 'summary' function for map objects.
