\name{print}
\alias{print}
\alias{print.summary.map}
\title{Print a Map Summary Object}
\description{Generate formatted output for objects of type 'summary.map'.}
\usage{
\S3method{print}{summary.map}(x,...)
}
\arguments{
  \item{x}{An object of type 'summary.map'.}
    \item{...}{Not used.}
}

\author{Lutz Hamel}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map(df,labels,xdim=15,ydim=10,train=)

## create a summary object and print it
s <- summary(m)
print(s)
}
