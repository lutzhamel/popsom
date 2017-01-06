\name{map.neuron}
\alias{map.neuron}
\title{Return a Map Neuron}
\description{Returns the contents of a neuron at (x,y) on the map as a vector.}
\usage{
map.neuron(map, x, y)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{x}{map x-coordinate of neuron.}
  \item{y}{map y-coordinate of neuron.}
}
\value{A vector representing the neuron.}

\author{Lutz Hamel}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display the neuron at position (1,1)
map.neuron(m,1,1)
}