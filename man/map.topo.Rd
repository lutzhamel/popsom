\name{map.topo}
\alias{map.topo}
\title{Estimated Topographical Accuracy}
\description{Evaluate the topological quality of a SOM using the estimated topographical accuracy.}
\usage{
map.topo(map,k=50,conf.int = 0.95,verb=FALSE,interval=TRUE)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{k}{number of samples to use in the computation of the estimated topographical accuracy (default=50)}
  \item{conf.int}{the confidence interval of the estimated topographical accuracy (default 95 percent).}
  \item{verb}{a switch controlling the structure of the output value (default=FALSE)}
  \item{interval}{a switch that controls whether the confidence interval is computed (default=TRUE)}

}  
\value{
	1) The value of the estimated topographical accuracy.
	2) The low value of the confidence interval and the high value of the confidence interval 'conf.int' if interval=TRUE.  
	3) If verb=TRUE then map.accuracy will return a vector with the accuracies of the individual k samples.
		}
\note{map.topo has been deprecated as of this release (version 4.3.0) of the popsom package for R. map.topo has also been flagged for removal upon release of the next version of the popsom package.

The use of map.convergence, also available in the popsom package, with the argument 'verb = True' is now encouraged to explicitly compute the estimated topographic accuracy of maps in future popsom package releases.
}
\references{
"SOM Quality Measures: A Statistical Approach," 
Lutz Hamel, WSOM16, 2016'. 
}
\author{Lutz Hamel}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display estimated topographical accuracy of the map
map.topo(m)

}
