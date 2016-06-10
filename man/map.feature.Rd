\name{map.feature}
\alias{map.feature}
\title{Generate Enhanced Unified Distance Matrix For Feature}
\description{Computes and displays the enhanced unified distance matrix for a feature of the training data.}
\usage{
map.feature(map, feature, explicit = FALSE, smoothing = 2)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{feature}{an integer as the index of the feature.}
  \item{explicit}{controls the shape of the connected components.}
  \item{smoothing}{controls the smoothing level of the umat (NULL, 0, >0).}
}
\note{This essentially constructs a horizontal slice of the UMAT.}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display the umat for the first feature of the map
map.feature(m, 1)
}