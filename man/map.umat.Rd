\name{map.umat}
\alias{map.umat}
\title{Generate Unified Distance Matrix For Map}
\description{Computes and displays the unified distance matrix for the map.} 
\usage{
map.umat(map)
}
\arguments{
  \item{map}{an object of type 'map'.}
}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display the umat for the map
map.umat(m)
}