### map-utils.R
# (c) University of Rhode Island
#               Lutz Hamel, Benjamin Ott, Greg Breard,
#               Robert Tatoian, Vishakh Gopu, Michael Eiger
#
# This file constitues a set of routines which are useful in constructing
# and evaluating self-organizing maps (SOMs).
# The main utilities available in this file are:
# map.build -------- constructs a map
# map.summary ------ compute a summary object
# map.convergence -- details of the convergence index of a map
# map.starburst ---- displays the starburst representation of the SOM model,
#                    the centers of starbursts are the centers of clusters
# map.fitted ------- returns a vector of labels assigned to the observations
# map.predict ------ returns classification labels for points in DF
# map.position ----- return the position of points on the map
# map.significance - graphically reports the significance of each feature with
#                    respect to the self-organizing map model
# map.marginal ----- displays a density plot of a training dataframe dimension
#                    overlayed with the neuron density for that same dimension or
#                    index.
#
### License
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
### Helpful docs
# medium.com/@ODSC/unsupervised-learning-evaluating-clusters-bd47eed175ce
# en.wikipedia.org/wiki/Total_sum_of_squares
###

# load libraries
loadNamespace("fields")
loadNamespace("graphics")
loadNamespace("ggplot2")
loadNamespace("hash")

### constructor ###

# map.build -- construct a SOM, returns an object of class 'map'
#
# parameters:
# - data - a dataframe where each row contains an unlabeled training instance
# - labels - a vector or dataframe with one label for each observation in data
# - xdim,ydim - the dimensions of the map
# - alpha - the learning rate, should be a positive non-zero real number
# - train - number of training iterations
# - normalize - normalize the input data by row
# - seed - a seed value for repeatablity of random initialization and selection
# - minimal - when true only the trained neuron's are returned.
# value:
# - an object of type 'map' -- see below
map.build <- function(data,
                      labels=NULL,
                      xdim=10,
                      ydim=5,
                      alpha=.3,
                      train=1000,
                      normalize=FALSE,
                      seed=NULL,
                      minimal=FALSE)
{
  if (alpha <= 0 || alpha > 1)
    stop("invalid value for alpha")

  if (xdim < 5 || ydim < 5)
    stop("map is too small.")

  if (!is.data.frame(data))
    stop("training data has to be a data frame")

  if (!all(sapply(data,is.numeric)))
    stop("only numeric data can be used for training")

  if (!is.null(labels) && !is.data.frame(labels))
    labels <- data.frame(labels)

  if (normalize)
    data <- map.normalize(data)

  if (!is.null(seed) && seed <= 0)
    stop("seed value has to be a positive integer value")

  if (!is.null(seed) && !test.integer(seed))
    stop("seed value has to be a positive integer value")

  if (!test.integer(train))
    stop("train value has to be a positive integer value")

  # train the neural network
  neurons <- vsom.f(data,
                    xdim=xdim,
                    ydim=ydim,
                    alpha=alpha,
                    train=train,
                    seed=seed)

  # make the neuron data a data frame
  neurons <- data.frame(neurons)
  names(neurons) <- names(data)

  # construct the map object
  map <- list(data=data,
              labels=labels,
              xdim=xdim,
              ydim=ydim,
              alpha=alpha,
              train=train,
              normalize=normalize,
              seed=seed,
              neurons=neurons)

  # add the class name
  if (minimal)
  {
    class(map) <- "map.minimal"
    return(map)
  }
  else
  {
    class(map) <- "map"
  }

  # NOTE: do not change the order of the following computations

  # add the heat map
  map$heat <- compute.heat(map)

  # list of indexes of the best matching neuron for each observation.
  # each index is an row index into the 'neuron' data frame.
  map$fitted.obs <- map.fitted.obs(map)

  # map of centroids - this is a map where each cell points
  # to the the location on the map where the corresponding centroid
  # is located.  centroids point to themselves.
  map$centroids <- compute.centroids(map)

  # a list of actual centroid locations on the map
  map$unique.centroids <- get.unique.centroids(map)

  # this is a map where locations of centroids have a label associated
  # with them. if unlabeled data then we invented labels for the centroids
  map$centroid.labels <- majority.labels(map)

  # label-to-centroid lookup table
  # Note: a label can be associated with multiple centroids
  map$label.to.centroid <- compute.label.to.centroid(map)

  # a vector of lists of observations per centroid indexed
  # by the centroid number from unique.centroids
  map$centroid.obs <- compute.centroid.obs(map)

  ### quality measures of the model ###

  # the convergence index of a map
  map$convergence <- map.convergence(map,verb=FALSE)

  # compute the average within cluster sum of squares (wcss)
  # this is the average distance variance within the clusters of the
  # underlying cluster model
  map$wcss <- compute.wcss(map)

  # compute the between cluster sum of squares (bcss)
  # this is the distance variance between the cluster centroids of the
  # underlying cluster model
  map$bcss <- compute.bcss(map)

  return(map)
}

# map.summary - compute a summary object
# parameters:
#   map - an object of type 'map'
#   verb -  a switch controlling the output
# value:
#   a summary object of type 'summary.map'
map.summary <- function(map, verb=TRUE)
{
  if (class(map) != "map")
    stop("first argument is not a map object.")

  value <- list()

  # training parameters
  header <- c("xdim",
              "ydim",
              "alpha",
              "train",
              "normalize",
              "seed",
              "instances",
              "columns")
  v <- c(map$xdim,
         map$ydim,
         map$alpha,
         map$train,
         if (map$normalize) "TRUE" else "FALSE",
         if (is.null(map$seed)) "NULL" else map$seed,
         nrow(map$data),
         ncol(map$data))
  df <- data.frame(t(v))
  names(df) <- header
  row.names(df) <- " "
  value$training.parameters <- df

  # quality assessments
  header <- c("convergence","separation","clusters")
  v <- c(map$convergence,
         1.0 - map$wcss/map$bcss,
         length(map$unique.centroids))
  df <- data.frame(t(v))
  names(df) <- header
  row.names(df) <- " "
  value$quality.assessments <- df

  class(value) <- "summary.map"

  if (verb)
  {
      cat("\n")
      cat("Training Parameters:\n")
      print(value$training.parameters)
      cat("\n")

      cat("Quality Assessments:\n")
      print(format(value$quality.assessments,digits=2))
      cat("\n")
  }
  else
  {
      value
  }
}

# map.starburst - compute and display the starburst representation of clusters
# parameters:
# - map is an object of type 'map'
map.starburst <- function(map)
{
  if (class(map) != "map")
    stop("first argument is not a map object.")

  plot.heat(map)
}

# map.significance - compute the relative significance of each feature and plot it
# parameters:
# - map is an object if type 'map'
# - graphics is a switch that controls whether a plot is generated or not
# - feature.labels is a switch to allow the plotting of feature names vs
#   feature indices
# return value:
# - a vector containing the significance for each feature
map.significance <- function(map,graphics=TRUE,feature.labels=TRUE)
{
  if (class(map) != "map")
    stop("first argument is not a map object.")

  data.df <- data.frame(map$data)
  nfeatures <- ncol(data.df)

  # Compute the variance of each feature on the map
  var.v <- array(data=1,dim=nfeatures)
  for (i in 1:nfeatures)
  {
    var.v[i] <- var(data.df[[i]])
  }

  # we use the variance of a feature as likelihood of
  # being an important feature, compute the Bayesian
  # probability of significance using uniform priors

  var.sum <- sum(var.v)
  prob.v <- var.v/var.sum

  # plot the significance
  if (graphics)
  {
    par.v <- map.graphics.set()

    y <- max(prob.v)
    plot.new()
    plot.window(xlim=c(1,nfeatures),ylim=c(0,y))
    box()

    title(xlab="Features",ylab="Significance")

    xticks <- seq(1,nfeatures,1)
    yticks <- seq(0,y,y/4)
    if (feature.labels)
      xlabels <- names(data.df)
    else
      xlabels <- seq(1,nfeatures,1)
    ylabels <- formatC(seq(0,y,y/4),digits=2)
    axis(1,at=xticks,labels=xlabels)
    axis(2,at=yticks,labels=ylabels)

    points(1:nfeatures,prob.v,type="h")

    map.graphics.reset(par.v)
    }
    else
    {
    prob.v
  }
}

# map.marginal - plot that shows the marginal probability distribution of the
#                neurons and data
# parameters:
# - map is an object of type 'map'
# - marginal is the name of a training data frame dimension or index
map.marginal <- function(map,marginal)
{
  # ensure that map is a 'map' object
  if (class(map) != "map")
    stop("first argument is not a map object.")

  # check if the second argument is of type character
  if (!typeof(marginal) == "character")
  {
    train <- data.frame(points = map$data[[marginal]])
    neurons <- data.frame(points = map$neurons[[marginal]])

    train$legend <- 'training data'
    neurons$legend <- 'neurons'

    hist <- rbind(train,neurons)
    ggplot(hist, aes(points, fill = legend)) +
         geom_density(alpha = 0.2) +
         xlab(names(map$data)[marginal])
  }
  else if (marginal %in% names(map$data))
  {
    train <- data.frame(points = map$data[names(map$data) == marginal])
    colnames(train) <- c("points")

    neurons <- data.frame(points = map$neurons[names(map$neurons) == marginal])
    colnames(neurons) <- c("points")

    train$legend <- 'training data'
    neurons$legend <- 'neurons'

    hist <- rbind(train,neurons)
    ggplot(hist, aes(points, fill = legend)) +
         geom_density(alpha = 0.2) +
         xlab(marginal)
  }
  else
  {
    stop("second argument is not a data dimension or index")
  }
}

# map.fitted -- returns a vector of labels assigned to the observations
# parameters:
# - map is an object of type 'map'
# value:
# - a vector of labels
map.fitted <- function(map)
{
  if (class(map) != "map")
    stop("first argument is not a map object.")

  nobs <- length(map$fitted.obs)
  labels <- c()
  for (i in 1:nobs)
  {
    nix <- map$fitted.obs[[i]]
    coord <- coordinate(map,nix)
    x <- coord$x
    y <- coord$y
    c.x <- map$centroids[[x,y]]$x
    c.y <- map$centroids[[x,y]]$y
    l <- map$centroid.labels[[c.x,c.y]]
    labels <- c(labels,l)
  }
  labels
}

# map.predict -- returns classification labels for points in DF
# parameters:
# - map -- map object
# - points  -- data frame of points to be classified
# value:
# - the label of the centroid x belongs to
map.predict<- function (map,points)
{
  # local function to do the actual prediction
  predict.point <- function (x)
  {
    if (!is.vector(x))
      stop("argument has to be a vector.")

    if (length(x) != ncol(map$data))
      stop("predict vector dimensionality is incompatible")

    if (map$normalize)
      x <- as.vector(map.normalize(x))

    # find best matching neuron
    nix <- best.match(map,x)

    # find corresponding centroid
    coord <- coordinate(map,nix)
    ix <- coord$x
    iy <- coord$y
    c.xy <- map$centroids[[ix,iy]]
    c.nix <- rowix(map,c.xy)
    label <- map$centroid.labels[[c.xy$x,c.xy$y]]
    c.ix <- find.centroidix(map,c.xy)

    # compute the confidence of the prediction
    # compute x to centroid distance
    vectors <- rbind(map$neurons[c.nix,],x)
    x.to.c.distance <- max(as.matrix(dist(vectors))[1,])

    # compute the max radius of cluster
    # NOTE: we are using the training data here NOT the neurons
    vectors <- map$neurons[c.nix,]
    for (i in 1:length(map$centroid.obs[[c.ix]]))
    {
      obs.ix <- map$centroid.obs[[c.ix]][i]
      vectors <- rbind(vectors,map$data[obs.ix,])
    }
    max.o.to.c.distance <- max(as.matrix(dist(vectors))[1,])
    # add a little bit of slack so we don't wind up with a 0 confidence value
    max.o.to.c.distance <- max.o.to.c.distance + 0.05*max.o.to.c.distance

    # compute confidence value
    conf <- 1.0 - (x.to.c.distance/max.o.to.c.distance)
    return (c(label,conf))
  }

  if (is.vector(points))
    points <- t(data.frame(points))

  m <- data.frame(t(apply(points,1,predict.point)))
  names(m) <- c("Label", "Confidence")
  m
}

# map.position-- return the position of points on the map
# parameters:
# - map -- map object
# - points   -- a data frame of points to be mapped
# value:
# - x-y coordinates of points in points
map.position <- function (map,points)
{
  # local function to positon a point on the map
  position.point <- function(x)
  {
    if (!is.vector(x))
      stop("argument has to be a vector.")

    if (length(x) != ncol(map$data))
      stop("vector dimensionality is incompatible")

    if (map$normalize)
      x <- as.vector(map.normalize(x))

    nix <- best.match(map,x)
    coord <- coordinate(map,nix)
    return (c(coord$x,coord$y))
  }

  if (is.vector(points))
    points <- t(data.frame(points))

  m <- data.frame(t(apply(points,1,position.point)))
  names(m) <- c("x-dim", "y-dim")
  m
}

# map.convergence - details of the convergence index of a map
#
# parameters:
# - map is an object if type 'map'
# - conf.int is the confidence interval of the quality assessment (default 95%)
# - k is the number of samples used for the estimated topographic accuracy
#   computation
# - verb if true reports the two convergence components separately, otherwise
#   it will report the linear combination of the two
# - ks is a switch, true for ks-test, false for standard var and means test
#
# - return value is the convergence index
map.convergence <- function(map,conf.int=.95,k=50,verb=TRUE,ks = TRUE)
{
    if (ks)
      embed <- map.embed.ks(map,conf.int,verb=FALSE)
    else
      embed <- map.embed.vm(map,conf.int,verb=FALSE)

    topo <- map.topo(map,k,conf.int,verb=FALSE,interval=FALSE)

    if (verb)
        return (list(embed=embed,topo=topo))
    else
        return (0.5*embed + 0.5*topo)
}

#############################################################################
############################# local functions ###############################
#############################################################################

# test.integer -- test to see if x is an integer value
test.integer <- function(x)
{
  all.equal(x, as.integer(x), check.attributes = FALSE)
}

# find.centroidix -- given a coordinate find the ix into the
#                    unique.centroids table
find.centroidix <- function (map,cd)
{
  for (i in 1:length(map$unique.centroids))
  {
    if (cd$x == map$unique.centroids[[i]]$x &&
        cd$y == map$unique.centroids[[i]]$y)
      {
        return (i)
      }
  }
  stop("coordinate not a centroid")
}

# compute.centroid.obs -- compute the observations that belong to each
#                         centroid.
compute.centroid.obs <- function (map)
{
  centroid.obs <- array(list(),dim=length(map$unique.centroids))

  for (cluster.ix in 1:length(map$unique.centroids))
  {
    c.nix <- rowix(map,map$unique.centroids[[cluster.ix]])
    for (i in 1:nrow(map$data))
    {
      # find the centroid of the current observation's
      # best matching neuron
      coord <- coordinate(map,map$fitted.obs[i])
      # centroid of cluster the neuron belongs to
      c.obj.nix <- rowix(map,map$centroids[[coord$x,coord$y]])
      # if observation centroid equal current centroid add to vectors
      if (c.obj.nix == c.nix)
      {
        centroid.obs[[cluster.ix]] <- append(centroid.obs[[cluster.ix]],i)
      }
    }
  }
  as.vector(centroid.obs)
}


# compute.wcss -- compute the average within cluster sum of squares
# see here:
# medium.com/@ODSC/unsupervised-learning-evaluating-clusters-bd47eed175ce
compute.wcss <- function (map)
{
  # for each cluster gather all the point vectors that belong to that
  # cluster into table 'vectors' making sure that the centroid vector
  # is always the first vector in the table.  Then compute the
  # sum square distances from the centroid to all the points.
  # when computing the average make sure that we ignore the centroid vector.
  clusters.ss <- c()
  for (cluster.ix in 1:length(map$unique.centroids))
  {
    c.nix <- rowix(map,map$unique.centroids[[cluster.ix]])
    vectors <- map$neurons[c.nix,]
    for (i in 1:length(map$centroid.obs[[cluster.ix]]))
    {
      obs.ix <- map$centroid.obs[[cluster.ix]][i]
      vectors <- rbind(vectors,map$data[obs.ix,])
    }
    distances <- as.matrix(dist(vectors))[1,]
    distances.sqd <- sapply(distances,function(x){x*x})
    c.ss <- sum(distances.sqd)/(length(distances.sqd)-1)
    clusters.ss <- c(clusters.ss,c.ss)
  }
  wcss <- sum(clusters.ss)/length(clusters.ss)
  wcss
}

# compute.bcss -- compute the average between cluster sum of squares
# see here:
# medium.com/@ODSC/unsupervised-learning-evaluating-clusters-bd47eed175ce
compute.bcss <- function (map)
{
  all.bc.ss <- c()

  # put all cluster vectors into one table
  c.nix <- rowix(map,map$unique.centroids[[1]])
  cluster.vectors <- map$neurons[c.nix,]
  if (length(map$unique.centroids) > 1)
  {
    for (cluster.ix in 2:length(map$unique.centroids))
    {
      c.nix <- rowix(map,map$unique.centroids[[cluster.ix]])
      c.vector <- map$neurons[c.nix,]
      cluster.vectors <- rbind(cluster.vectors,c.vector)
    }
  }

  # put each cluster vector at the beginning of the table in turn
  # and compute the distances - row 1 will have our results
  # NOTE: at every iteration one of the cluster vectors will
  # appear twice in the vector table. we have to adjust for that
  # when computing the average.
  for (cluster.ix in 1:length(map$unique.centroids))
  {
    c.nix <- rowix(map,map$unique.centroids[[cluster.ix]])
    c.vector <- map$neurons[c.nix,]
    compute.vectors <- rbind(c.vector,cluster.vectors)
    bc.distances <- as.matrix(dist(compute.vectors))[1,]
    bc.distances.sqd <- sapply(bc.distances,function(x){x*x})
    # Note: cluster.ix appears twice
    bc.ss <- sum(bc.distances.sqd)/(length(bc.distances.sqd)-2)
    all.bc.ss <- c(all.bc.ss,bc.ss)
  }

  # basic sanity check
  stopifnot(length(all.bc.ss)==length(map$unique.centroids))
  bcss <- sum(all.bc.ss)/length(all.bc.ss)
  bcss
}

# compute.label.to.centroid -- compute a label to centroid lookup table
# The returned value is a table of indexes into the unique centroids table
compute.label.to.centroid <- function (map)
{
  conv <- hash()

  for (i in 1:length(map$unique.centroids))
  {
    x <- map$unique.centroids[[i]]$x
    y <- map$unique.centroids[[i]]$y
    l <- map$centroid.labels[[x,y]]
    if (is.null(conv[[l]]))
    {
      conv[[l]] <- list(i)
    }
    else
    {
      conv[[l]] <- append(conv[[l]],i)
    }
  }
  conv
}

# for each observation i, visual has an entry for
# the best matching neuron
map.fitted.obs <- function(map)
{
  fitted.obs <- c()
  for (i in 1:nrow(map$data))
  {
      b <- best.match(map,map$data[i,])
      fitted.obs <- c(fitted.obs,b)
  }

  fitted.obs
}

# map.topo - measure the topographic accuracy of the map using sampling
#
# parameters:
# - map is an object if type 'map'
# - k is the number of samples used for the accuracy computation
# - conf.int is the confidence interval of the accuracy test (default 95%)
# - verb is switch that governs the return value, false: single accuracy value
#   is returned, true: a vector of individual feature accuracies is returned.
# - interval is a switch that controls whether the confidence interval
#   is computed.
#
# - return value is the estimated topographic accuracy
map.topo <- function(map,k=50,conf.int=.95,verb=FALSE,interval=TRUE)
{
  if (class(map) != "map")
      stop("first argument is not a map object.")

  # data.df is a matrix that contains the training data
  data.df <- as.matrix(map$data)

  # sample map$data
  # k samples unless k > nrow(data.df), then
  k <- min(k,nrow(data.df))
  data.sample.ix <- sample(1:nrow(data.df),size=k,replace=FALSE)

  # compute the sum topographic accuracy - the accuracy of a single sample
  # is 1 if the best matching unit is a neighbor otherwise it is 0
  acc.v <- c()
  for (i in 1:k)
  {
      acc.v <- c(acc.v,
                 accuracy(map,
                          data.df[data.sample.ix[i],],
                          data.sample.ix[i]))
  }

  # compute the confidence interval values using the bootstrap
  if (interval)
      bval <- bootstrap(map,conf.int,data.df,k,acc.v)

  # the sum topographic accuracy is scaled by the number of samples - estimated
  # topographic accuracy
  if (verb)
  {
      acc.v
  }
  else
  {
      val <- sum(acc.v)/k
      if (interval)
          list(val=val,lo=bval$lo,hi=bval$hi)
      else
          val

  }
}

# map.embed using variance and mean tests
map.embed.vm <- function(map,conf.int=.95,verb=FALSE)
{

    if (class(map) != "map")
        stop("first argument is not a map object.")

    # map.df is a dataframe that contains the neurons
    map.df <- data.frame(map$neurons)

    # data.df is a dataframe that contain the training data
    # note: map$data is what the 'som' package returns
    data.df <- data.frame(map$data)

    # do the F-test on a pair of datasets: code vectors/training data
    vl <- df.var.test(map.df,data.df,conf=conf.int)

    # do the t-test on a pair of datasets: code vectors/training data
    ml <- df.mean.test(map.df,data.df,conf=conf.int)

    # compute the variance captured by the map -- but only if the
    # means have converged as well.
    nfeatures <- ncol(map.df)
    prob.v <- map.significance(map,graphics=FALSE)
    var.sum <- 0

    for (i in 1:nfeatures)
    {
        if (vl$conf.int.lo[i] <= 1.0
            && vl$conf.int.hi[i] >= 1.0
            &&  ml$conf.int.lo[i] <= 0.0
            && ml$conf.int.hi[i] >= 0.0)
        {
            var.sum <- var.sum + prob.v[i]
        }
        else
        {
            # not converged - zero out the probability
            prob.v[i] <- 0
        }
    }

    # return the variance captured by converged features
    if (verb)
      prob.v
    else
      var.sum
}

# map.embed using the kolgomorov-smirnov test
map.embed.ks <- function(map,conf.int=.95,verb=FALSE)
{
  if (class(map) != "map")
  {
      stop("first argument is not a map object.")
  }

  # map.df is a dataframe that contains the neurons
  map.df <- data.frame(map$neurons)

  # data.df is a dataframe that contain the training data
  # note: map$data is what the 'som' package returns
  data.df <- data.frame(map$data)

  nfeatures <- ncol(map.df)

  # use the Kolmogorov-Smirnov Test to test whether the neurons and
  # training data appear to come from the same distribution
  ks.vector <- NULL
  for(i in 1:nfeatures){
      # suppress the warning about ties.
      ks.vector[[i]] <- suppressWarnings(ks.test(map.df[[i]], data.df[[i]]))
  }

  prob.v <- map.significance(map,graphics=FALSE)
  var.sum <- 0

  # compute the variance captured by the map
  for (i in 1:nfeatures)
  {
      # the second entry contains the p-value
      if (ks.vector[[i]][[2]] > (1 - conf.int)) {
          var.sum <- var.sum + prob.v[i]
      } else {
          # not converged - zero out the probability
          prob.v[i] <- 0
      }
  }

  # return the variance captured by converged features
  if (verb)
    prob.v
  else
    var.sum
}

# map.normalize -- based on the som:normalize function but preserved names
map.normalize <- function (x)
{
  if (is.vector(x))
  {
    return (scale(x))
  }
  else if (is.data.frame(x))
  {
    df <- data.frame(t(apply(x,1,scale)))
    names(df) <- names(x)
    return (df)
  }
  else
  {
    stop("'x' is not a vector or dataframe.\n")
  }
}


# bootstrap -- compute the topographic accuracies for the given
#              confidence interval
bootstrap <- function(map,conf.int,data.df,k,sample.acc.v)
{
  ix <- as.integer(100 - conf.int*100)
  bn <- 200

  bootstrap.acc.v <- c(sum(sample.acc.v)/k)

  for (i in 2:bn)
  {
      bs.v <- sample(1:k,size=k,replace=TRUE)
      a <- sum(sample.acc.v[bs.v])/k
      bootstrap.acc.v <- c(bootstrap.acc.v,a)
  }

  bootstrap.acc.sort.v <- sort(bootstrap.acc.v)

  lo.val <- bootstrap.acc.sort.v[ix]
  hi.val <- bootstrap.acc.sort.v[bn-ix]

  list(lo=lo.val,hi=hi.val)
}

# best.match -- given observation obs, return the best matching neuron
best.match <- function(map,obs,full=FALSE)
{
  # NOTE: replicate obs so that there are nr rows of obs
  obs.m <- matrix(as.numeric(obs),
                  nrow(map$neurons),
                  ncol(map$neurons),
                  byrow=TRUE)
  diff <- map$neurons - obs.m
  squ <- diff * diff
  s <- rowSums(squ)
  d <- sqrt(s)
  o <- order(d)

  if (full)
    o
  else
    o[1]
}

# accuracy -- the topographic accuracy of a single sample is 1 is the best
#             matching unit and the second best matching unit are are neighbors
#             otherwise it is 0
accuracy <- function(map,sample,data.ix)
{
    # compute the euclidean distances of the sample from the neurons
    # and find the 2 best matching units for the sample

    o <- best.match(map,sample,full=TRUE)
    best.ix <- o[1]
    second.best.ix <- o[2]

    # sanity check
    coord <- coordinate(map,best.ix)
    coord.x <- coord$x
    coord.y <- coord$y

    map.ix <- map$fitted.obs[data.ix]
    coord <- coordinate(map,map.ix)
    map.x <- coord$x
    map.y <- coord$y

    if (coord.x != map.x || coord.y != map.y || best.ix != map.ix)
    {
        cat("best.ix: ",best.ix," map.rix: ",map.ix,"\n")
        stop("problems with coordinates")
    }

    # determine if the best and second best are neighbors on the map
    best.xy <- coordinate(map,best.ix)
    second.best.xy <- coordinate(map,second.best.ix)
    diff.map <- c(best.xy$x,best.xy$y) - c(second.best.xy$x,second.best.xy$y)
    diff.map.sq <- diff.map * diff.map
    sum.map <- sum(diff.map.sq)
    dist.map <- sqrt(sum.map)

    # it is a neighbor if the distance on the map
    # between the bmu and 2bmu is less than 2
    if (dist.map < 2)
      1
    else
      0
}

# coord -- constructor for a 'coord' object
coord <- function (x=-1,y=-1)
{
  l <- list(x=x,y=y)
  class(l) <- "coord"
  l
}

# coordinate -- convert from a row index to a map xy-coordinate
coordinate <- function(map,rowix)
{
    x <- (rowix-1) %% map$xdim + 1
    y <- (rowix-1) %/% map$xdim + 1
    coord(x,y)
}

# rowix -- convert from a map xy-coordinate to a row index
rowix <- function(map,cd)
{
    if (class(cd) != "coord")
      stop("expected a coord object")

    rix <- cd$x + (cd$y-1)*map$xdim
    rix
}

# map.graphics.set -- set the graphics environment for our map utilities
#                     the return value is the original graphics param vector
map.graphics.set <- function()
{
  par.v <- par()
  par(ps=6)
  par.v
}

# map.graphics.reset -- reset the graphics environment to the original state
# parameter - a vector containing the settings for the original state
map.graphics.reset <- function(par.vector)
{
  par(ps=par.vector$ps)
}

# plot.heat - plot a heat map based on a 'map', this plot also contains the
#             connected components of the map based on the landscape of the
#             heat map
plot.heat <- function(map)
{
  x <- map$xdim
  y <- map$ydim
  centroids <- map$centroids

  ### need to make sure the map doesn't have a dimension of 1
  if (x <= 1 || y <= 1)
  {
    stop("map dimensions too small")
  }

  ### bin the heat values into 100 bins used for the 100 heat colors
  heat.v <- as.vector(map$heat)
  heat.v <- cut(heat.v,breaks=100,labels=FALSE)
  heat <- array(data=heat.v,dim=c(x,y))
  colors<- heat.colors(100)

  ### set up the graphics window
  par.v <- map.graphics.set()
  plot.new()
  plot.window(xlim=c(0,x),ylim=c(0,y))
  box()

  title(xlab="x",ylab="y")

  xticks <- seq(0.5,x-0.5,1)
  yticks <- seq(0.5,y-0.5,1)
  xlabels <- seq(1,x,1)
  ylabels <- seq(1,y,1)
  axis(1,at=xticks,labels=xlabels)
  axis(3,at=xticks,labels=xlabels)
  axis(2,at=yticks,labels=ylabels)
  axis(4,at=yticks,labels=ylabels)

  ### plot the neurons as heat squares on the map
  # TODO: vectorize this - rect can operate on vectors of coordinates and values
  for (ix in 1:x)
  {
    for (iy in 1:y)
    {
      rect(ix-1,iy-1,ix,iy,col=colors[100 - heat[ix,iy] + 1],border=NA)
    }
  }

  # connect each neuron to its centroid
  for(ix in 1:x)
  {
    for (iy in 1:y)
    {
      cx <- centroids[[ix,iy]]$x
      cy <- centroids[[ix,iy]]$y
      points(c(ix,cx)-.5,c(iy,cy)-.5,type="l",col="grey")
    }
  }

  # put majority labels on the centroids
  # Note: if labels were not given then the function majority.labels
  # will compute numerical labels to attach to the centroids.
  centroid.labels <- majority.labels(map)

  for(ix in 1:x)
  {
    for (iy in 1:y)
    {
      lab <- centroid.labels[[ix,iy]]
      if (lab != none.label)
      {
        text(ix-.5,iy-.5,labels=lab)
      }
    }
  }

  map.graphics.reset(par.v)
}


#compute.centroids -- compute the centroid for each point on the map
# parameters:
# - map is an object if type 'map'
# return value:
# - a map of the same dimension of the input map where each cell points
#   to the centroid of this cell.
compute.centroids <- function(map)
{
  heat <- map$heat
  xdim <- map$xdim
  ydim <- map$ydim
  max.val <- max(heat)
  centroids <- array(data=list(coord()),dim=c(xdim,ydim))

  ########################################################################
  ### local recursive function to find the centroid of a point on the map
  compute.centroid <- function(ix,iy)
  {
    # first we check if the current position is already associated
    # with a centroid.  if so, simply return the coordinates
    # of that centroid
    if ((centroids[[ix,iy]])$x > -1 && (centroids[[ix,iy]])$y > -1)
    {
      centroids[[ix,iy]]
    }

    # try to find a smaller value in the immediate neighborhood
    # make our current position the square with the minimum value.
    # if a minimum value other than our own current value cannot be
    # found then we are at a minimum.
    #
    # search the neighborhood; three different cases: inner element,
    # corner element, side element
    # TODO: there has to be a better way!

    min.val <- heat[ix,iy]
    min.x <- ix
    min.y <- iy

    # (ix,iy) is an inner map element
    if (ix > 1 && ix < xdim && iy > 1 && iy < ydim)
    {
      if (heat[ix-1,iy-1] < min.val)
      {
        min.val <- heat[ix-1,iy-1]
        min.x <- ix-1
        min.y <- iy-1
      }
      if (heat[ix,iy-1] < min.val)
      {
        min.val <- heat[ix,iy-1]
        min.x <- ix
        min.y <- iy-1
      }
      if (heat[ix+1,iy-1] < min.val)
      {
        min.val <- heat[ix+1,iy-1]
        min.x <- ix+1
        min.y <- iy-1
      }
      if (heat[ix+1,iy] < min.val)
            {
        min.val <- heat[ix+1,iy]
        min.x <- ix+1
        min.y <- iy
      }
      if (heat[ix+1,iy+1] < min.val)
            {
        min.val <- heat[ix+1,iy+1]
        min.x <- ix+1
        min.y <- iy+1
      }
      if (heat[ix,iy+1] < min.val)
      {
        min.val <- heat[ix,iy+1]
        min.x <- ix
        min.y <- iy+1
      }
      if (heat[ix-1,iy+1] < min.val)
      {
        min.val <- heat[ix-1,iy+1]
        min.x <- ix-1
        min.y <- iy+1
      }
      if (heat[ix-1,iy] < min.val)
      {
        min.val <- heat[ix-1,iy]
        min.x <- ix-1
        min.y <- iy
      }
    }

    # (ix,iy) is bottom left corner
    else if (ix == 1 && iy == 1)
    {
      if (heat[ix+1,iy] < min.val)
            {
        min.val <- heat[ix+1,iy]
        min.x <- ix+1
        min.y <- iy
      }
      if (heat[ix+1,iy+1] < min.val)
            {
        min.val <- heat[ix+1,iy+1]
        min.x <- ix+1
        min.y <- iy+1
      }
      if (heat[ix,iy+1] < min.val)
            {
        min.val <- heat[ix,iy+1]
        min.x <- ix
        min.y <- iy+1
      }
    }

    # (ix,iy) is bottom right corner
    else if (ix == xdim && iy == 1)
    {
      if (heat[ix,iy+1] < min.val)
      {
        min.val <- heat[ix,iy+1]
        min.x <- ix
        min.y <- iy+1
      }
      if (heat[ix-1,iy+1] < min.val)
      {
        min.val <- heat[ix-1,iy+1]
        min.x <- ix-1
        min.y <- iy+1
      }
      if (heat[ix-1,iy] < min.val)
      {
        min.val <- heat[ix-1,iy]
        min.x <- ix-1
        min.y <- iy
      }
    }

    # (ix,iy) is top right corner
    else if (ix == xdim && iy == ydim)
    {
      if (heat[ix-1,iy-1] < min.val)
      {
        min.val <- heat[ix-1,iy-1]
        min.x <- ix-1
        min.y <- iy-1
      }
      if (heat[ix,iy-1] < min.val)
      {
        min.val <- heat[ix,iy-1]
        min.x <- ix
        min.y <- iy-1
      }
      if (heat[ix-1,iy] < min.val)
      {
        min.val <- heat[ix-1,iy]
        min.x <- ix-1
        min.y <- iy
      }
    }

    # (ix,iy) is top left corner
    else if (ix == 1 && iy == ydim)
    {
      if (heat[ix,iy-1] < min.val)
      {
        min.val <- heat[ix,iy-1]
        min.x <- ix
        min.y <- iy-1
      }
      if (heat[ix+1,iy-1] < min.val)
      {
        min.val <- heat[ix+1,iy-1]
        min.x <- ix+1
        min.y <- iy-1
      }
      if (heat[ix+1,iy] < min.val)
      {
        min.val <- heat[ix+1,iy]
        min.x <- ix+1
        min.y <- iy
      }
    }

    # (ix,iy) is a left side element
    else if (ix == 1  && iy > 1 && iy < ydim)
    {
      if (heat[ix,iy-1] < min.val)
      {
        min.val <- heat[ix,iy-1]
        min.x <- ix
        min.y <- iy-1
      }
      if (heat[ix+1,iy-1] < min.val)
      {
        min.val <- heat[ix+1,iy-1]
        min.x <- ix+1
        min.y <- iy-1
      }
      if (heat[ix+1,iy] < min.val)
      {
        min.val <- heat[ix+1,iy]
        min.x <- ix+1
        min.y <- iy
      }
      if (heat[ix+1,iy+1] < min.val)
      {
        min.val <- heat[ix+1,iy+1]
        min.x <- ix+1
        min.y <- iy+1
      }
      if (heat[ix,iy+1] < min.val)
      {
        min.val <- heat[ix,iy+1]
        min.x <- ix
        min.y <- iy+1
      }
    }

    # (ix,iy) is a bottom side element
    else if (ix > 1 && ix < xdim && iy == 1 )
    {
      if (heat[ix+1,iy] < min.val)
      {
        min.val <- heat[ix+1,iy]
        min.x <- ix+1
        min.y <- iy
      }
      if (heat[ix+1,iy+1] < min.val)
      {
        min.val <- heat[ix+1,iy+1]
        min.x <- ix+1
        min.y <- iy+1
      }
      if (heat[ix,iy+1] < min.val)
      {
        min.val <- heat[ix,iy+1]
        min.x <- ix
        min.y <- iy+1
      }
      if (heat[ix-1,iy+1] < min.val)
      {
        min.val <- heat[ix-1,iy+1]
        min.x <- ix-1
        min.y <- iy+1
      }
      if (heat[ix-1,iy] < min.val)
      {
        min.val <- heat[ix-1,iy]
        min.x <- ix-1
        min.y <- iy
      }
    }

    # (ix,iy) is a right side element
    else if (ix == xdim && iy > 1 && iy < ydim)
    {
      if (heat[ix-1,iy-1] < min.val)
      {
        min.val <- heat[ix-1,iy-1]
        min.x <- ix-1
        min.y <- iy-1
      }
      if (heat[ix,iy-1] < min.val)
      {
        min.val <- heat[ix,iy-1]
        min.x <- ix
        min.y <- iy-1
      }
      if (heat[ix,iy+1] < min.val)
      {
        min.val <- heat[ix,iy+1]
        min.x <- ix
        min.y <- iy+1
      }
      if (heat[ix-1,iy+1] < min.val)
      {
        min.val <- heat[ix-1,iy+1]
        min.x <- ix-1
        min.y <- iy+1
      }
      if (heat[ix-1,iy] < min.val)
      {
        min.val <- heat[ix-1,iy]
        min.x <- ix-1
        min.y <- iy
      }
    }

    # (ix,iy) is a top side element
    else if (ix > 1 && ix < xdim && iy == ydim)
    {
      if (heat[ix-1,iy-1] < min.val)
      {
        min.val <- heat[ix-1,iy-1]
        min.x <- ix-1
        min.y <- iy-1
      }
      if (heat[ix,iy-1] < min.val)
      {
        min.val <- heat[ix,iy-1]
        min.x <- ix
        min.y <- iy-1
      }
      if (heat[ix+1,iy-1] < min.val)
      {
        min.val <- heat[ix+1,iy-1]
        min.x <- ix+1
        min.y <- iy-1
      }
      if (heat[ix+1,iy] < min.val)
      {
        min.val <- heat[ix+1,iy]
        min.x <- ix+1
        min.y <- iy
      }
      if (heat[ix-1,iy] < min.val)
      {
        min.val <- heat[ix-1,iy]
        min.x <- ix-1
        min.y <- iy
      }
    }

    # if successful
    # move to the square with the smaller value and
    # call compute.centroid on this new square
    if (min.x != ix || min.y != iy)
    {
      # Note: returns a list of an x and y coordinate
      compute.centroid(min.x,min.y)
    }
    #else
    # we have found a minimum -- this is our centroid.
    else
    {
      coord(ix,iy)
    }
  } # end function compute.centroid
  ###########################################################################

  ### iterate over the map and find the centroid for each element
  for (i in 1:xdim)
  {
    for (j in 1:ydim)
    {
      centroids[[i,j]] <- compute.centroid(i,j)
    }
  }
  centroids
}

# compute.heat -- compute a heat value map representation of the given
#                 distance matrix
# parameters:
# - map is an object if type 'map'
# return value:
# - a matrix with the same x-y dims as the original map containing the heat
compute.heat <- function(map)
{
  d <- as.matrix(dist(data.frame(map$neurons)))
  x <- map$xdim
  y <- map$ydim
  heat <- array(data=0,dim=c(x,y))

  if (x == 1 || y == 1)
    stop("heat map cannot be computed for a map of dimension 1")

  # local function as a shorthand for rowix
  xl <- function(ix,iy)
  {
    #cat("converting (",ix,",",iy,") to row", ix + (iy-1) *xdim,"\n")
    #ix + (iy-1) * x
    rowix(map,coord(ix,iy))
  }

  # check if the map is larger than 2 x 2 (otherwise it is only corners)
  if (x > 2 && y > 2)
  {
    # iterate over the inner nodes and compute their umat values
    for (ix in 2:(x-1))
    {
      for (iy in 2:(y-1))
      {
        sum <-
             d[xl(ix,iy),xl(ix-1,iy-1)] +
             d[xl(ix,iy),xl(ix,iy-1)] +
             d[xl(ix,iy),xl(ix+1,iy-1)] +
             d[xl(ix,iy),xl(ix+1,iy)] +
             d[xl(ix,iy),xl(ix+1,iy+1)] +
             d[xl(ix,iy),xl(ix,iy+1)] +
             d[xl(ix,iy),xl(ix-1,iy+1)] +
             d[xl(ix,iy),xl(ix-1,iy)]
        heat[ix,iy] <- sum/8
      }
    }

    # iterate over bottom x axis
    for (ix in 2:(x-1))
    {
      iy <- 1
      sum <-
           d[xl(ix,iy),xl(ix+1,iy)] +
           d[xl(ix,iy),xl(ix+1,iy+1)] +
           d[xl(ix,iy),xl(ix,iy+1)] +
           d[xl(ix,iy),xl(ix-1,iy+1)] +
           d[xl(ix,iy),xl(ix-1,iy)]
      heat[ix,iy] <- sum/5
    }

    # iterate over top x axis
    for (ix in 2:(x-1))
    {
      iy <- y
      sum <-
           d[xl(ix,iy),xl(ix-1,iy-1)] +
           d[xl(ix,iy),xl(ix,iy-1)] +
           d[xl(ix,iy),xl(ix+1,iy-1)] +
           d[xl(ix,iy),xl(ix+1,iy)] +
           d[xl(ix,iy),xl(ix-1,iy)]
      heat[ix,iy] <- sum/5
    }

    # iterate over the left y-axis
    for (iy in 2:(y-1))
    {
      ix <- 1
      sum <-
           d[xl(ix,iy),xl(ix,iy-1)] +
           d[xl(ix,iy),xl(ix+1,iy-1)] +
           d[xl(ix,iy),xl(ix+1,iy)] +
           d[xl(ix,iy),xl(ix+1,iy+1)] +
           d[xl(ix,iy),xl(ix,iy+1)]
      heat[ix,iy] <- sum/5
    }

    # iterate over the right y-axis
    for (iy in 2:(y-1))
    {
      ix <- x
      sum <-
           d[xl(ix,iy),xl(ix-1,iy-1)] +
           d[xl(ix,iy),xl(ix,iy-1)] +
           d[xl(ix,iy),xl(ix,iy+1)] +
           d[xl(ix,iy),xl(ix-1,iy+1)] +
           d[xl(ix,iy),xl(ix-1,iy)]
      heat[ix,iy] <- sum/5
    }
  } # end if

  # compute umat values for corners
  if (x >= 2 && y >= 2)
    {
    # bottom left corner
    ix <- 1
    iy <- 1
    sum <-
        d[xl(ix,iy),xl(ix+1,iy)] +
        d[xl(ix,iy),xl(ix+1,iy+1)] +
        d[xl(ix,iy),xl(ix,iy+1)]
    heat[ix,iy] <- sum/3

    # bottom right corner
    ix <- x
    iy <- 1
    sum <-
         d[xl(ix,iy),xl(ix,iy+1)] +
         d[xl(ix,iy),xl(ix-1,iy+1)] +
         d[xl(ix,iy),xl(ix-1,iy)]
    heat[ix,iy] <- sum/3

    # top left corner
    ix <- 1
    iy <- y
    sum <-
        d[xl(ix,iy),xl(ix,iy-1)] +
        d[xl(ix,iy),xl(ix+1,iy-1)] +
        d[xl(ix,iy),xl(ix+1,iy)]
    heat[ix,iy] <- sum/3

    # top right corner
    ix <- x
    iy <- y
    sum <-
        d[xl(ix,iy),xl(ix-1,iy-1)] +
        d[xl(ix,iy),xl(ix,iy-1)] +
        d[xl(ix,iy),xl(ix-1,iy)]
    heat[ix,iy] <- sum/3
  } # end if

  # smooth the heat map
  xcoords <- c()
  ycoords <- c()
  for (i in 1:y)
  {
    for (j in 1:x)
    {
      ycoords <- c(ycoords, i)
      xcoords <- c(xcoords, j)
    }
  }
  xycoords <- data.frame(xcoords,ycoords)
  heat <- smooth.2d(as.vector(heat),
                    x=as.matrix(xycoords),
                    nrow=x,ncol=y,
                    surface=FALSE,
                    theta=2)
  heat
}

# df.var.test -- a function that applies the F-test testing the ratio
#                  of the variances of the two data frames
# parameters:
# - df1,df2 - data frames with the same number of columns
# - conf - confidence level for the F-test (default .95)
df.var.test <- function(df1,df2,conf = .95)
{
  if (length(df1) != length(df2))
        stop("cannot compare variances of data frames")

  # init our working arrays
  var.ratio.v <- array(data=1,dim=length(df1))
  var.confintlo.v <- array(data=1,dim=length(df1))
  var.confinthi.v <- array(data=1,dim=length(df1))

  # compute the F-test on each feature in our populations
  for (i in 1:length(df1))
    {
    t <- var.test(df1[[i]],df2[[i]],conf.level=conf)
    var.ratio.v[i] <- t$estimate
    #cat("Feature",i,"confidence interval =",t$conf.int,"\n")
    var.confintlo.v[i] <- t$conf.int[1]
    var.confinthi.v[i] <- t$conf.int[2]
  }

  # return a list with the ratios and conf intervals for each feature
  list(ratio=var.ratio.v,
       conf.int.lo=var.confintlo.v,
       conf.int.hi=var.confinthi.v)
}

# df.mean.test -- a function that applies the t-test testing the difference
#                   of the means of the two data frames
# parameters:
# - df1,df2 - data frames with the same number of columns
# - conf - confidence level for the t-test (default .95)
df.mean.test <- function(df1,df2,conf = .95)
{
  if (ncol(df1) != ncol(df2))
        stop("cannot compare means of data frames")

  # init our working arrays
  mean.diff.v <- array(data=1,dim=ncol(df1))
  mean.confintlo.v <- array(data=1,dim=ncol(df1))
  mean.confinthi.v <- array(data=1,dim=ncol(df1))

  # compute the F-test on each feature in our populations
  for (i in 1:ncol(df1)) {
    t <- t.test(x=df1[[i]],y=df2[[i]],conf.level=conf)
    mean.diff.v[i] <- t$estimate[1] - t$estimate[2]
    #cat("Feature",i,"confidence interval =",t$conf.int,"\n")
    mean.confintlo.v[i] <- t$conf.int[1]
    mean.confinthi.v[i] <- t$conf.int[2]
  }

  # return a list with the mean differences and conf intervals for each feature
  list(diff=mean.diff.v,
       conf.int.lo=mean.confintlo.v,
       conf.int.hi=mean.confinthi.v)
}


# vsom.f - vectorized and optimized version of the stochastic
# SOM training algorithm written in Fortran90
vsom.f <- function(data,xdim,ydim,alpha,train,seed)
{
    ### some constants
    dr <- nrow(data)
    dc <- ncol(data)
    nr <- xdim*ydim
    nc <- dc # dim of data and neurons is the same

    ### build and initialize the matrix holding the neurons
    if (!is.null(seed))
    {
        set.seed(seed)
    }
    else
    {
        # send a -1 to Fortran function to indicate no seed present
        seed <- -1
    }
    cells <- nr * nc        # no. of neurons times number of data dimensions
    v <- runif(cells,-1,1)  # vector with small init values for all neurons
    # NOTE: each row represents a neuron, each column represents a dimension.
    neurons <- matrix(v,nrow=nr,ncol=nc)  # rearrange the vector as matrix

    result <- .Fortran("vsom",
                       as.single(neurons),
                       as.single(data.matrix(data)),
                       as.integer(dr),
                       as.integer(dc),
                       as.integer(xdim),
                       as.integer(ydim),
                       as.single(alpha),
                       as.integer(train),
                       as.integer(seed),
                       PACKAGE="popsom")

    # unpack the structure and list in result[1]
    v <- result[1]
    # rearrange the result vector as matrix
    neurons <- matrix(v[[1]],nrow=nr,ncol=nc,byrow=FALSE)

    neurons
}

# get.unique.centroids -- a list of unique centroid locations on the map
get.unique.centroids <- function(map)
{
  # get the dimensions of the map
  centroids <- map$centroids
  xdim <- map$xdim
  ydim <- map$ydim
  cd.list <- c()
  for(ix in 1:xdim)
  {
    for(iy in 1:ydim)
    {
      c.xy <- centroids[[ix, iy]]
      b <- sapply(cd.list, function (x) {x$x == c.xy$x && x$y == c.xy$y})
      if (!any(b))
      {
        cd.list <- c(cd.list,list(c.xy))
      }
    }
  }
  as.vector(cd.list)
}


# majority.labels -- return a map where the positions of the centroids
# has the majority label of the appropriate cluster attached to them.

none.label <- "<None>"

majority.labels <- function(map)
{
  if (is.null(map$labels))
  {
    # no labels given, make up some numerical labels
    return (numerical.labels(map))
  }

  x <- map$xdim
  y <- map$ydim
  centroids <- map$centroids
  nobs <- nrow(map$data)
  centroid.labels <- array(data=list(),dim=c(x,y))
  majority.labels <- array(data=none.label,dim=c(x,y))

  # gather the labels from the clusters and record them
  # at the centroid position.
  for(i in 1:nobs)
  {
   lab <- as.character(map$labels[i,1])
   nix <- map$fitted.obs[i]
   c <- coordinate(map,nix)
   ix <- c$x
   iy <- c$y
   cx <- centroids[[ix,iy]]$x
   cy <- centroids[[ix,iy]]$y
   centroid.labels[[cx,cy]] <- append(centroid.labels[[cx,cy]],lab)
  }

  # attach majority labels to centroids
  for (ix in 1:x)
  {
   for (iy in 1:y)
   {
     label.v <- centroid.labels[[ix,iy]]
     if (length(label.v)!=0)
     {
       majority <- data.frame(sort(table(label.v),decreasing=TRUE))
       if (nrow(majority) == 1) # only one label
       {
         # just copy a label from the label vector
         majority.labels[[ix,iy]] <- label.v[1]
       }
       else
       {
         majority.labels[[ix,iy]] <- levels(majority[1,1])[1]
       }
     }
   }
  }
  majority.labels
}

# numerical.labels -- create labels for centroids
numerical.labels <- function(map)
{
  label_cnt <- 1
  centroids <- map$centroids
  unique.centroids <- map$unique.centroids
  centroid.labels <- array(data=none.label,dim=c(map$xdim,map$ydim))

  # set our labels at the centroid locations
  for (i in 1:length(unique.centroids))
  {
    # create a label
    label <- paste("centroid",label_cnt,sep=" ")
    label_cnt <- label_cnt+1
    ix <- unique.centroids[[i]]$x
    iy <- unique.centroids[[i]]$y
    #cat("coord",ix,iy,label,"\n")
    centroid.labels[[ix,iy]] <- label
  }
  centroid.labels
}

########################## research stuff ###########################

# compute.nwcss -- compute the average within cluster sum of squares
#                  of neuron clusters
compute.nwcss <- function (map)
{
  # for each cluster gather all the point vectors that belong to that
  # cluster into table 'vectors' making sure that the centroid vector
  # is always the first vector in the table.  Then compute the
  # sum square distances from the centroid to all the points.
  # when computing the average make sure that we ignore the centroid vector.
  clusters.ss <- c()
  for (cluster.ix in 1:length(map$unique.centroids))
  {
    c.nix <- rowix(map,map$unique.centroids[[cluster.ix]])
    vectors <- map$neurons[c.nix,]
    for (i in 1:length(map$centroid.obs[[cluster.ix]]))
    {
      obs.ix <- map$centroid.obs[[cluster.ix]][i]
      obs.nix <- map$fitted[[obs.ix]]
      obs.coord <- coordinate(map,obs.nix)
      centroid.coord <- map$centroids[[obs.coord$x,obs.coord$y]]
      centroid.nix <- rowix(map,centroid.coord)
      if (centroid.nix == c.nix){
        vectors <- rbind(vectors,map$neurons[obs.nix,])
      }
    }
    distances <- as.matrix(dist(vectors))[1,]
    distances.sqd <- sapply(distances,function(x){x*x})
    c.ss <- sum(distances.sqd)/(length(distances.sqd)-1)
    clusters.ss <- c(clusters.ss,c.ss)
  }
  wcss <- sum(clusters.ss)/length(clusters.ss)
  wcss
}

# avg.homogeneity -- given labels another way to ascertain quality of the map
avg.homogeneity <- function(map)
{
  if (is.null(map$labels))
  {
    stop("you need to attach labels to the map")
  }

  # need to make sure the map doesn't have a dimension of 1
  if (map$xdim <= 1 || map$ydim <= 1)
  {
    stop("map dimensions too small")
  }

  x <- map$xdim
  y <- map$ydim
  centroids <- map$centroids
  nobs <- nrow(map$data)
  centroid.labels <- array(data=list(),dim=c(x,y))

  #attach labels to centroids
  # count the labels in each map cell
  for(i in 1:nobs)
  {
   lab <- as.character(map$labels[i,1])
   nix <- map$fitted.obs[i]
   c <- coordinate(map,nix)
   ix <- c$x
   iy <- c$y
   cx <- centroids[[ix,iy]]$x
   cy <- centroids[[ix,iy]]$y
   centroid.labels[[cx,cy]] <- append(centroid.labels[[cx,cy]],lab)
  }

  # compute average homogeneity of the map: h = (1/nobs)*sum_c majority.label_c
  sum.majority <- 0
  n.centroids <- 0

  for (ix in 1:x)
  {
   for (iy in 1:y)
   {
     label.v <- centroid.labels[[ix,iy]]
     if (length(label.v)!=0)
     {
       n.centroids <- n.centroids + 1
       majority <- data.frame(sort(table(label.v),decreasing=TRUE))

       if (nrow(majority) == 1) # only one label
       {
         m.val <- length(label.v)
       }
       else
       {
         m.val <- majority[1,2]
       }
       sum.majority <- sum.majority + m.val
     }
   }
  }
  list(homog=sum.majority/nobs, nclust=n.centroids)
}
