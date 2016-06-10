\name{map.accuracy}
\alias{map.accuracy}
\title{Estimated Topographical Accuracy}
\description{Evaluate the topological quality of a SOM using the estimated topographical accuracy.}
\usage{
map.accuracy(map, k=50, conf.int = 0.95, verb=FALSE)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{k}{number of samples to use in the computation of the estimated topographical accuracy (default=50)}
  \item{conf.int}{the confidence interval of the estimated topographical accuracy (default 95 percent).}
  \item{verb}{a switch controlling the structure of the output value (default=FALSE)}
}  
\value{
	A list with three components: 1) the value of the estimated topographical accuracy 2) the low value
	of the confidence interval 'conf.int' 3) the high value of the confidence interval.  If verb=TRUE then
	map.accuracy will return a vector with the accuracies of the individual k samples.
		}
\references{
"SOM Quality Measures: A Statistical Approach," 
Lutz Hamel, in preparation, 2015. 
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
map.accuracy(m)

}