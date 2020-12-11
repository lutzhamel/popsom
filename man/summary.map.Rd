\name{summary}
\alias{summary}
\alias{summary.map}
\title{Summary Object for Map}
\description{Generate a summary object for 'map' objects.}
\usage{
\S3method{summary}{map}(object,...)
}
\arguments{
  \item{object}{An object of type 'map'.}
    \item{...}{Not used.}
}

\value{
An object of type 'summary.map' which contains two structures:
\describe{
\item{training.parameters}{A dataframe containing the parameters the map was trained with.}
\item{quality.assessments}{A dataframe containing the quality assessments of the map. In particular, it
    contains the 'convergence'  of the map which is a linear combination of variance capture  and
    topographic fidelity of the map.  Furthermore, it contains the 'spread' of the clusters indicating
    the separation of the clusters.  This is computed by the formula,
    \preformatted{ 1 - wcss/bcss}
    In general, a value close to 1 means well separated clusters.
}
}
}

\author{Lutz Hamel}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map(df,labels,xdim=15,ydim=10,train=10000)

## compute a summary object and display it
s <- summary(m)
s
}
