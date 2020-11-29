\name{map}
\alias{map}
\title{Build Map}
\usage{
map(data,labels=NULL,xdim=10,ydim=5,alpha=0.3,train=1000,normalize=TRUE)
}
\arguments{
  \item{data}{A dataframe where each row contains an unlabeled training instance.}
  \item{labels}{A vector or dataframe with one label for each observation in data.}
  \item{xdim}{The x-dimension of the map. (default=10)}
  \item{ydim}{The y-dimension of the map. (default=5)}
  \item{alpha}{The learning rate, should be a positive non-zero real number. (default=0.3)}
  \item{train}{The number of training iterations. (default=1000)}
  \item{normalize}{Boolean switch indicating whether or not to normalize the data. (default=TRUE)}
}
\description{The function 'map' constructs an object of type 'map'. The object contains two models: (1) A self-organizing map model expressed through its trained neurons and its quality of fit can be ascertained by the 'convergence' (see below). (2) A cluster model expressed by the discovered centroids and its quality can be ascertained by the 'within cluster sum of squares' and the 'between cluster sum of squares' (see below).}
\value{An object of type 'map'.  The object has the following member fields:
\describe{
\item{data}{Data frame contining the possibly normalized training data.}
\item{labels}{Vector of labels, one for each observation in data or NULL if no labels were given.}
\item{xdim}{The x dimension of the neuron map.}
\item{ydim}{The y dimension of the neuron map.}
\item{alpha}{The given learning rate for the neural network.}
\item{train}{The training iterations applied to the neural network.}
\item{neurons}{A list of neurons for the network. The dimensionality of this
data frame is the same as the training data.  The following two formulas
come in handy when working with the neural data.  The first set of computations
provide the x and y coordinate on the map of the neuron in row 'rowix' of the 'neurons' data frame,
\preformatted{
    x <- (rowix-1) %% map$xdim + 1
    y <- (rowix-1) %/% map$xdim + 1
}
The second formula computes the row of the neuron in position (x,y) on the map,
\preformatted{
    rowix <- x + (y-1)*map$xdim
}}
\item{heat}{This is the representation of the map which is the basis for the 'starburst' plot.}
\item{fitted.obs}{List of indexes of the best matching neuron for each observation.  Each index is an row index into the 'neuron' data frame.}
\item{centroids}{This is a data frame of (x,y) locations where each cell points
  to the the (x,y) location on the map where the corresponding centroid
  is located.  Centroids point to themselves.}
\item{unique.centroids}{A vector of actual centroid (x,y)-locations on the map. Hint: to compute the number of clusters on the map take the length of this vector.}
\item{centroid.labels}{A data frame where the (x,y)-locations of actual centroids have a label associated with them. All other locations are NULL. If the training data is unlabeled  then we invented a label for each of the centroids.}
\item{label.to.centroid}{A label-to-centroid lookup table (hash). A lookup in this table will return a list of indexes into the 'unique.centroids' table. Note: a label can be associated with multiple centroids.}
\item{centroid.obs}{A vector of lists of observations per centroid indexed
by the centroid number from 'unique.centroids'.  The observations on the list are row numbers of the 'data' data frame.}
\item{convergence}{A quality measure of how well the map fit the training data.}
\item{wcss}{The average 'within cluster sum of squares'. This is the average
distance variance within the clusters of the underlying cluster model.}
\item{bcss}{The 'between cluster sum of squares'.  This is the distance 
variance between the cluster centroids of the underlying cluster model.}
}
}

\note{If your training data is unlabeled popsom will automatically
generate a label for each of the centroids it discovers.}
\author{Lutz Hamel, Benjamin Ott, Gregory Breard}
\examples{
data(iris)

df <- subset(iris,select=-Species)
labels <- subset(iris,select=Species)

m <- map(df, labels, xdim=15, ydim=10, train=1000)
}
