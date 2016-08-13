\name{map.starburst}
\alias{map.starburst}
\title{Generate Starburst For Map}
\description{Computes and displays the starburst representation of clusters on a map.}
\usage{
map.starburst(map, explicit = FALSE, smoothing = 2)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{explicit}{controls the shape of the connected components.}
  \item{smoothing}{controls the smoothing level of the map display (NULL, 0, >0).}
}

\references{
"Improved Interpretability of the Unified Distance Matrix with Connected Components,"
Lutz Hamel and Chris W. Brown. Proceeding of the 7th International Conference on 
Data Mining (DMIN'11), July 18-21, 2011, Las Vegas Nevada, USA, ISBN: 1-60132-168-6, pp338-343, 
CSREA Press, 2011.
}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display the starburst for the map
map.starburst(m)
}