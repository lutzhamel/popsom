\name{fitted}
\alias{fitted}
\alias{fitted.map}
\title{Fit Observations}
\description{Computes a vector of labels assigned to each of the observations in the training data through the constructed cluster model.  If the training data is unlabeled then generated labels are used.}
\usage{
\S3method{fitted}{map}(map)
}
\arguments{
  \item{map}{An object of type 'map'.}
}
\value{A vector of predicted labels, one for each observations in the training data.}


\author{Lutz Hamel}
\examples{
data(iris)

df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

m <- map(df,labels,xdim=15,ydim=10,train=10000)

fitted(m)
}
