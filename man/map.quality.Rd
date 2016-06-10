\name{map.quality}
\alias{map.quality}
\title{SOM Quality Assessment}
\description{Evaluate the quality of a SOM using the embedding and estimated topographical accuracy.}
\usage{
map.quality(map,conf.int=.95,k=50)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{conf.int}{is the confidence interval of the quality assessment (default 95\%)}
  \item{k}{number of samples to use in the computation of the estimated topographical accuracy (default=50)}
}  
\value{
	A pair of values: 1) embedding accuracy 2) estimated topographic accuracy
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

## map quality
map.quality(m)

}