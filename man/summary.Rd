\name{summary}
\alias{summary}
\title{Summary Object for Map}
\description{Generate a summary object for 'map' objects.}
\usage{
summary(map,verb=TRUE)
}
\arguments{
  \item{map}{An object of type 'map'.}
  \item{verb}{A switch controlling the output.}
}

\value{
An object of type 'summary.map' which contains two structures:
\describe{
\item{training.parameters}{A dataframe containing the parameters the map was trained with.}
\item{quality.assessments}{A dataframe containing the quality assessments of the map. In particular, it
    contains the 'convergence'  of the map which is a linear combination of variance capture  and
    topographic fidelity of the map. A value close to 1 means a converged map (for more details
    see the reference below).  Furthermore, it contains the 'separation' of the clusters.
    This is computed by the formula,
    \preformatted{ 1 - wcss/bcss}
    In general, a value close to 1 means well separated clusters.}
}
If 'verb' is TRUE the summar.map object will be formatted and printed to the screen, otherwise
it will be returned as a data structure.
}

\author{Lutz Hamel}
\references{
Self-Organizing Map Convergence, Robert Tatoian and Lutz Hamel. Proceedings of the 2016 International Conference on Data Mining (DMIN'16), pp92-98, July 25-28, 2016, Las Vegas, Nevada, USA, ISBN: 1-60132-431-6, CSREA Press.
}
\examples{
data(iris)

## set data frame and labels
df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

## build a map
m <- map(df,labels,xdim=15,ydim=10,train=10000)

## compute a summary object and display it
s <- summary(m,verb=FALSE)
s
}
