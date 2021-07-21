\name{predict}
\alias{predict}
\title{Compute Classification Labels for Given Points}
\description{Compute classification labels for points in a given data frame using the underlying clustering model.  If the training data is unlabeled then generated labels are used.}
\usage{
predict(map,points)
}
\arguments{
  \item{map}{An object of type 'map'.}
  \item{points}{A data frame of points to be classified.}
}
\value{A data frame with classification results.  The data frame has two columns:
\describe{
\item{Label}{The assigned label to the observation at the same row in the 'points' data frame.}
\item{Confidence}{A confidence value assigned to the label prediction.}
}
}

\author{Lutz Hamel}
\examples{
data(iris)

df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

m <- map(df,labels,xdim=15,ydim=10,train=10000)

predict(m,df)
}
