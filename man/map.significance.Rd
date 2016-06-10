\name{map.significance}
\alias{map.significance}
\title{Compute Significance Of Features}
\description{Computes the relative significance of each feature and plots it.}
\usage{
map.significance(map, graphics = TRUE, feature.labels = TRUE)
}
\arguments{
  \item{map}{an object of type 'map'.}
  \item{graphics}{a switch that controls whether a plot is generated or not.}
  \item{feature.labels}{a switch to allow the plotting of feature names vs feature indices.}
}
\value{if graphics=FALSE a vector containing the significance for each feature is returned.}
\note{We use a Bayesian approach to compute the relative significance of features based on 
variance.}

\references{
"Bayesian Probability Approach to Feature Significance for Infrared Spectra of Bacteria,"
Lutz Hamel, Chris W. Brown, Applied Spectroscopy, Volume 66, Number 1, 2012.
}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map.build(df, labels, xdim=15, ydim=10, train=1000)

## show the relative feature significance for each feature
data.frame(names(df),map.significance(m,graphics=FALSE))

## display the relative feature significance graphically
map.significance(m)
}