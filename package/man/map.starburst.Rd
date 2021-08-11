\name{map.starburst}
\alias{map.starburst}
\title{Generate Starburst For Map}

\usage{
map.starburst(map)
}
\arguments{
  \item{map}{An object of type 'map'}
}

\description{
Generate a starburst representation of the clusters on the heat map for the
self-organizing map model.
}

\references{
Improved Interpretability of the Unified Distance Matrix with Connected Components,
Lutz Hamel and Chris W. Brown. Proceeding of the 7th International Conference on
Data Mining (DMIN'11), July 18-21, 2011, Las Vegas Nevada, USA, ISBN: 1-60132-168-6, pp338-343,
CSREA Press, 2011.
}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard, Robert Tatoian, Vishakh Gopu}
\examples{
data(iris)

df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

m <- map.build(df,labels,xdim=15,ydim=10,train=100000)

map.starburst(m)
}
