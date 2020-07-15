\name{map.projection}
\alias{map.projection}
\title{Map Projection}
\description{Prints the association of labels with map elements.}
\usage{
map.projection(map)
}
\arguments{
  \item{map}{an object of type 'map'.}
}
\value{a dataframe containing the projection onto the map for each training observation.}
\note{map.projection has been deprecated as of this release (version 4.3.0) of the popsom package for R. map.projection has also been flagged for removal upon release of the next version of the popsom package.
}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display the label association for the map
map.projection(m)
}
