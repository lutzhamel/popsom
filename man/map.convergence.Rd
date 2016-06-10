\name{map.convergence}
\alias{map.convergence}
\title{Evaluate Map Convergence}
\description{Evaluates the convergence of a map using a two-sample statistical test.}
\usage{
map.convergence(map, conf.int = 0.95, verb=FALSE)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{conf.int}{the confidence interval of the convergence test (default 95 percent).}
  \item{verb}{a switch controlling the structure of the output value (default=FALSE)}
}  
\value{The convergence index of the map.
	If the switch verb=TRUE then a vector of the individual feature convergence indices is returned.
	}
\note{The convergence index is the variance of the trainig data captured by the map. 
Maps with low convergence index are typically not trustworthy. 
However, the precise cut-off depends on the noise level in your training data.
}
\references{
"A Population Based Convergence Criterion for Self-Organizing Maps," 
Lutz Hamel and Benjamin Ott. Proceeding of the 2012 International Conference on Data Mining (DMIN'12), 
pp98-104, July 16-19, 2012, Las Vegas Nevada, USA. 
}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## display the convergence index of the map
map.convergence(m)

## display the convergence indices of the individual features
data.frame(names(df),map.convergence(m,verb=TRUE))
}