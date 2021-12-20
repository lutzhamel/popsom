\name{map.convergence}
\alias{map.convergence}
\title{SOM Quality Assessment}
\description{Evaluate the quality of a SOM using embedding accuracy and estimated topographical accuracy.}
\usage{
map.convergence(map,conf.int=.95,k=50,verb=TRUE,ks=TRUE)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{conf.int}{is the confidence interval of the quality assessment.}
  \item{k}{number of samples to use in the computation of the estimated topographical accuracy.}
  \item{verb}{if true reports the two convergence components separately, otherwise it will report a linear combination of the two indices.}
  \item{ks}{if true uses the Kolmogorov-Smirnov convergence test otherwise a convergence test based on variance and means is performed.}
}
\value{
	A single value or a pair of values: 1) embedding accuracy 2) estimated topographic accuracy.
  The structure of the return value depends on the 'verb' switch.
		}
\references{
"SOM Quality Measures: A Statistical Approach,"
Lutz Hamel, WSOM16, 2016.
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
map.convergence(m)

}
