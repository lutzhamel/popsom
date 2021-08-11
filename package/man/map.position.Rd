\name{map.position}
\alias{map.position}
\title{Compute Map Positions for Given Points}
\description{Compute the (x,y)-positions of points on the map.}
\usage{
map.position(map,points)
}
\arguments{
  \item{map}{An object of type 'map'.}
  \item{points}{A data frame of points to be mapped.}
}
\value{A data frame with (x,y)-positions.  The data frame has two columns:
\describe{
\item{x-dim}{The x-position of the corresponding point in the 'points' data frame.}
\item{y-dim}{The y-position of the corresponding point in the 'points' data frame.}
}
}

\author{Lutz Hamel}
\examples{
data(iris)

df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

m <- map.build(df,labels,xdim=15,ydim=10,train=10000)

map.position(m,df)
}
