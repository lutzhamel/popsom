\name{map.significance}
\alias{map.significance}
\title{Compute Significance of Features}
\description{Computes the relative significance of each feature and plots it.}
\usage{
map.significance(map,graphics=TRUE,feature.labels=TRUE)
}
\arguments{
  \item{map}{An object of type 'map'.}
  \item{graphics}{A switch that controls whether a plot is generated or not.}
  \item{feature.labels}{A switch to allow the plotting of feature names vs feature indices.}
}
\value{If graphics=FALSE a vector containing the significance for each feature is returned.}

\note{We use a Bayesian approach to compute the relative significance of features based on
variance.}

\references{
Bayesian Probability Approach to Feature Significance for Infrared Spectra of Bacteria,
Lutz Hamel, Chris W. Brown, Applied Spectroscopy, Volume 66, Number 1, 2012.
}
\author{Lutz Hamel}
\examples{
data(iris)

df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

m <- map.build(df,labels,xdim=15,ydim=10,train=10000)

## show the relative feature significance for each feature
df <- data.frame(names(df),map.significance(m,graphics=FALSE))
names(df) <- c("Features","Significance")
df

## display the relative feature significance graphically
map.significance(m)
}
