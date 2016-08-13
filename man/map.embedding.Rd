\name{map.embedding}
\alias{map.embedding}
\title{Compute the Map Embedding Accuracy}
\description{Evaluates how well a map models the underlying training data distribution.}
\usage{
map.embedding(map, conf.int = 0.95, verb=FALSE)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{conf.int}{the confidence interval of the embedding test (default 95 percent).}
  \item{verb}{a switch controlling the structure of the output value (default=FALSE)}
}  
\value{The embedding accuracy of the map.
	If the switch verb=TRUE then a vector of the individual feature embedding accuracies are returned.
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

## display the embedding accuracy of the map
map.embedding(m)

## display the embedding accuracies of the individual features
data.frame(names(df),map.embedding(m,verb=TRUE))
}