### map-utils.R
# version 5.0.0
# (c) 2009-2020 Lutz Hamel, Benjamin Ott, Greg Breard, University of Rhode Island
#               with Robert Tatoian, Vishakh Gopu, and Michael Eiger
#
# This file constitues a set of routines which are useful in constructing
# and evaluating self-organizing maps (SOMs).
# The main utilities available in this file are:
# map.build --------- constructs a map
# map.convergence --- reports the map convergence index
# map.significance -- graphically reports the significance of each feature with
#                     respect to the self-organizing map model
# map.starburst ----- displays the starburst representation of the SOM model, the centers of
#                     starbursts are the centers of clusters
# map.marginal ------ displays a density plot of a training dataframe dimension overlayed
#                      with the neuron density for that same dimension or index.
#
### License
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
###

# load libraries
require(som)
require(class)
require(fields)
require(graphics)
require(ggplot2)

### map.build -- construct a SOM, returns an object of class 'map'
# parameters:
# - data - a dataframe where each row contains an unlabeled training instance
# - labels - a vector or dataframe with one label for each observation in data
# - xdim,ydim - the dimensions of the map
# - alpha - the learning rate, should be a positive non-zero real number
# - train - number of training iterations
# - algorithm - selection switch
# - normalize - normalize the input data by row
# retuns:
# - an object of type 'map' -- see below

# Hint: if your training data does not have any labels you can construct a
#       simple label dataframe as follows: labels <- data.frame(1:nrow(training.data))

# NOTE: default algorithm: "vsom" also available: "som", "experimental", "batchsom"

map.build <- function(data,
                      labels=NULL,
                      xdim=10,
                      ydim=5,
                      alpha=.3,
                      train=1000,
                      algorithm="vsom",
                      normalize=TRUE)
{
  # check if the dims are reasonable
	if (xdim < 3 || ydim < 3)
		stop("map.build: map is too small.")

  if (normalize)
    data <- map.normalize(data)

  # pmatch: som == 1, vsom == 2, experimental == 3, batchsom == 4
  algorithms = c("som","vsom","experimental","batchsom")

  # train the map - returns a list of neurons
  if (pmatch(algorithm,algorithms,nomatch=0) == 1) # som
  {
      # compute the initial neighborhood radius
      r <- sqrt(xdim^2 + ydim^2)

      m <- som(data,
               xdim=xdim,
               ydim=ydim,
               init="random",
               alpha=c(alpha,alpha),
               alphaType="linear",
               neigh="bubble",
               topol="rect",
               radius=c(r,r),
               rlen=c(1,train))

      # the 'som' package does something really annoying with attributes
      # for the neuron matrix, we get rid of that by casting the neurons
      # as a new matrix
      neurons <- matrix(m$code,xdim*ydim,ncol(data))
  }
  else if (pmatch(algorithm,algorithms,nomatch=0) == 2) # vsom
  {
      neurons <- vsom.f(data,
                        xdim=xdim,
                        ydim=ydim,
                        alpha=alpha,
                        train=train)
  }
  else if (pmatch(algorithm,algorithms,nomatch=0) == 3) # experimental
  {
      neurons <- vsom.r(data,
                        xdim=xdim,
                        ydim=ydim,
                        alpha=alpha,
                        train=train)
  }
  else if (pmatch(algorithm,algorithms,nomatch=0) == 4) # batchsom
  {
      # compute the initial neighborhood radius
      r <- sqrt(xdim^2 + ydim^2)

      m <- batchsom.private(data.matrix(data),
                            grid=somgrid(xdim,ydim,"rectangular"),
                            min.radius=1,
                            max.radius=r,
                            train=train,
                            "random",
                            "bubble")

      # extract the neurons
      neurons <- matrix(m$codes,xdim*ydim,ncol(data))
  }
  else
  {
      stop("map.build only supports 'som','vsom','experimental',and 'batchsom'")
  }

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
              algorithm=algorithm,
              normalize=normalize,
              neurons=neurons)

  # add the class name
  class(map) <- "map"

  # add the heat map
  map$heat <- compute.heat(map)

  # ix's of best matching neuron for each obs
  map$visual <- map.visual(map)

  return(map)
}

### map.convergence - the convergence index of a map
#
# parameters:
# - map is an object if type 'map'
# - conf.int is the confidence interval of the quality assessment (default 95%)
# - k is the number of samples used for the estimated topographic accuracy computation
# - verb if true reports the two convergence components separately, otherwise it will
#        report the linear combination of the two
# - ks is a switch, true for ks-test, false for standard var and means test
#
# - return value is the convergence index

map.convergence <- function(map,conf.int=.95,k=50,verb=FALSE,ks = TRUE)
{
    if (ks)
      embed <- map.embed.ks(map,conf.int,verb=FALSE)
    else
      embed <- map.embed.vm(map,conf.int,verb=FALSE)

    topo <- map.topo.private(map,k,conf.int,verb=FALSE,interval=FALSE)

    if (verb)
        return (list(embed=embed,topo=topo))
    else
        return (0.5*embed + 0.5*topo)
}


### map.starburst - compute and display the starburst representation of clusters
# parameters:
# - map is an object if type 'map'
# - explicit controls the shape of the connected components

map.starburst <- function(map,explicit=FALSE)
{

	if (class(map) != "map")
		stop("map.starburst: first argument is not a map object.")

	plot.heat(map,explicit=explicit,comp=TRUE)
}

### map.significance - compute the relative significance of each feature and plot it
# parameters:
# - map is an object if type 'map'
# - graphics is a switch that controls whether a plot is generated or not
# - feature.labels is a switch to allow the plotting of feature names vs feature indices
# return value:
# - a vector containing the significance for each feature

map.significance <- function(map,graphics=TRUE,feature.labels=TRUE)
{
	if (class(map) != "map")
		stop("map.significance: first argument is not a map object.")

	data.df <- data.frame(map$data)
	nfeatures <- ncol(data.df)

	# Compute the variance of each feature on the map
	var.v <- array(data=1,dim=nfeatures)
	for (i in 1:nfeatures)
  {
		var.v[i] <- var(data.df[[i]]);
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



### map.marginal - plot that shows the marginal probability distribution of the neurons and data
# parameters:
# - map is an object of type 'map'
# - marginal is the name of a training data frame dimension or index

map.marginal <- function(map,marginal)
{
  # ensure that map is a 'map' object
  if (class(map) != "map")
    stop("map.marginal: first argument is not a map object.")

  # check if the second argument is of type character
  if (!typeof(marginal) == "character")
  {
    train <- data.frame(points = map$data[[marginal]])
    neurons <- data.frame(points = map$neurons[[marginal]])

    train$legend <- 'training data'
    neurons$legend <- 'neurons'

    hist <- rbind(train,neurons)
    ggplot(hist, aes(points, fill = legend)) + geom_density(alpha = 0.2) + xlab(names(map$data)[marginal])

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
    ggplot(hist, aes(points, fill = legend)) + geom_density(alpha = 0.2) + xlab(marginal)

  }
  else
  {
    stop("map.marginal: second argument is not the name of a training data frame dimension or index")
  }

}

############################### local functions #################################

### map.neuron - returns the contents of a neuron at (x,y) on the map as a vector
# parameters:
#   map - the neuron map
#   x - map x-coordinate of neuron
#   y - map y-coordinate of neuron
# return value:
#   a vector representing the neuron

map.neuron <- function(map,x,y)
{
    if (class(map) != "map")
        stop("map.neuron: first argument is not a map object.")

    ix <- rowix(map,x,y)
    map$neurons[ix,]
}

# for each observation i, visual has an entry for
# the best matching neuron

map.visual <- function(map)
{
  visual <- c()
  for (i in 1:nrow(map$data))
  {
      b <- best.match(map,map$data[i,])
      visual <- c(visual,b)
  }

  visual
}

map.topo.private <- function(map,k=50,conf.int=.95,verb=FALSE,interval=TRUE)
{


    if (class(map) != "map")
        stop("map.topo: first argument is not a map object.")

    # data.df is a matrix that contains the training data
    data.df <- as.matrix(map$data)

    # sample map$data
    # TODO: think of something clever here rather than just aborting.
    if (k > nrow(data.df))
        stop("map.topo: sample larger than training data.")

    data.sample.ix <- sample(1:nrow(data.df),size=k,replace=FALSE)

    # compute the sum topographic accuracy - the accuracy of a single sample
    # is 1 if the best matching unit is a neighbor otherwise it is 0
    acc.v <- c()
    for (i in 1:k)
    {
        acc.v <- c(acc.v,accuracy(map,data.df[data.sample.ix[i],],data.sample.ix[i]))
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
        stop("map.embed: first argument is not a map object.")

    # map.df is a dataframe that contains the neurons
    map.df <- data.frame(map$neurons)

    # data.df is a dataframe that contain the training data
    # note: map$data is what the 'som' package returns
    data.df <- data.frame(map$data)

    # do the F-test on a pair of datasets: code vectors/training data
    vl <- df.var.test(map.df,data.df,conf=conf.int)

    # do the t-test on a pair of datasets: code vectors/training data
    ml <- df.mean.test(map.df,data.df,conf=conf.int)

    # compute the variance captured by the map -- but only if the means have converged as well.
    nfeatures <- ncol(map.df)
    prob.v <- map.significance(map,graphics=FALSE)
    var.sum <- 0

    for (i in 1:nfeatures)
    {
        #cat("Feature",i,"variance:\t",vl$ratio[i],"\t(",vl$conf.int.lo[i],"-",vl$conf.int.hi[i],")\n")
        #cat("Feature",i,"mean:\t",ml$diff[i],"\t(",ml$conf.int.lo[i],"-",ml$conf.int.hi[i],")\n")
        if (vl$conf.int.lo[i] <= 1.0 && vl$conf.int.hi[i] >= 1.0 &&  ml$conf.int.lo[i] <= 0.0 && ml$conf.int.hi[i] >= 0.0) {
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

# map.embed using the kolgomorov-smirnov test

map.embed.ks <- function(map,conf.int=.95,verb=FALSE)
{
  if (class(map) != "map")
  {
      stop("map.embed: first argument is not a map object.")
  }

  # map.df is a dataframe that contains the neurons
  map.df <- data.frame(map$neurons)

  # data.df is a dataframe that contain the training data
  # note: map$data is what the 'som' package returns
  data.df <- data.frame(map$data)

  nfeatures <- ncol(map.df)

  # use the Kolmogorov-Smirnov Test to test whether the neurons and training data appear
  # to come from the same distribution
  ks.vector <- NULL
  for(i in 1:nfeatures){
      # Note rpt - I needed to use suppress warnings to suppress the warning about ties.
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

map.normalize <- function (x, byrow = TRUE)
{
    if (is.data.frame(x))
    {
        if (byrow)
            df <- data.frame(t(apply(x, 1, scale)))
        else
            df <- data.frame(apply(x, 2, scale))

        names(df) <- names(x)
        df
    }
    else
        stop("map.normalize: 'x' is not a dataframe.\n")
}


# bootstrap -- compute the topographic accuracies for the given confidence interval

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
  obs.m <- matrix(as.numeric(obs),nrow(map$neurons),ncol(map$neurons),byrow=TRUE)
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

# accuracy -- the topographic accuracy of a single sample is 1 is the best matching unit
#             and the second best matching unit are are neighbors otherwise it is 0

accuracy <- function(map,sample,data.ix)
{
    # compute the euclidean distances of the sample from the neurons
    # and find the 2 best matching units for the sample

    o <- best.match(map,sample,full=TRUE)
    best.ix <- o[1]
    second.best.ix <- o[2]

    # sanity check
    coord <- coordinate(map,best.ix)
    coord.x <- coord[1]
    coord.y <- coord[2]

    map.ix <- map$visual[data.ix]
    coord <- coordinate(map,map.ix)
    map.x <- coord[1]
    map.y <- coord[2]

    if (coord.x != map.x || coord.y != map.y || best.ix != map.ix)
    {
        cat("best.ix: ",best.ix," map.rix: ",map.ix,"\n")
        stop("accuracy: problems with coordinates")
    }

    # determine if the best and second best are neighbors on the map
    best.xy <- coordinate(map,best.ix)
    second.best.xy <- coordinate(map,second.best.ix)
    diff.map <- best.xy - second.best.xy
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

# coordinate -- convert from a row index to a map xy-coordinate

coordinate <- function(map,rowix)
{
    x <- (rowix-1) %% map$xdim + 1
    y <- (rowix-1) %/% map$xdim + 1
    c(x,y)
}

#rowix -- convert from a map xy-coordinate to a row index

rowix <- function(map,x,y)
{
    rix <- x + (y-1)*map$xdim
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

### plot.heat - plot a heat map based on a 'map', this plot also contains the connected
#               components of the map based on the landscape of the heat map
# parameters:
# - map is an object if type 'map'
# - heat is a 2D heat map of the map returned by 'map'
# - labels is a vector with labels of the original training data set
# - explicit controls the shape of the connected components
# - comp controls whether we plot the connected components on the heat map

plot.heat <- function(map,explicit=FALSE,comp=TRUE)
{
	x <- map$xdim
	y <- map$ydim
	nobs <- nrow(map$data)
	count <- array(data=0,dim=c(x,y))

	### need to make sure the map doesn't have a dimension of 1
	if (x <= 1 || y <= 1)
  {
    stop("plot.heat: map dimensions too small")
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

	### put the connected component lines on the map
	if (comp)
  {
	  # find the centroid for each neuron on the map
	  centroids <- compute.centroids(map,explicit)
    # connect each neuron to its centroid
		for(ix in 1:x)
    {
			for (iy in 1:y)
      {
				cx <- centroids$centroid.x[ix,iy]
				cy <- centroids$centroid.y[ix,iy]
				points(c(ix,cx)-.5,c(iy,cy)-.5,type="l",col="grey")
			}
		}

    # put majority labels on the centroids
    if (!is.null(map$labels))
    {
      centroid.labels <- majority.labels(map)
      for(ix in 1:x)
      {
        for (iy in 1:y)
        {
          lab.l <- centroid.labels[[ix,iy]]
          if (length(lab.l)!=0)
          {
            text(ix-.5,iy-.5,labels=lab.l[1])
          }
        }
      }
    }
	}
	map.graphics.reset(par.v)
}


### compute.centroids -- compute the centroid for each point on the map
# parameters:
# - map is an object if type 'map'
# - explicit controls the shape of the connected component
# return value:
# - a list containing the matrices with the same x-y dims as the original map containing the centroid x-y coordinates

compute.centroids <- function(map,explicit=FALSE)
{
  heat <- map$heat
	xdim <- map$xdim
	ydim <- map$ydim
	centroid.x <- array(data=-1,dim=c(xdim,ydim))
	centroid.y <- array(data=-1,dim=c(xdim,ydim))
	max.val <- max(heat)

  ########################################################################
  ### local recursive function to find the centroid of a point on the map
	compute.centroid <- function(ix,iy)
  {
  	# first we check if the current position is already associated
  	# with a centroid.  if so, simply return the coordinates
  	# of that centroid
  	if (centroid.x[ix,iy] > -1 && centroid.y[ix,iy] > -1)
    {
  		list(bestx=centroid.x[ix,iy],besty=centroid.y[ix,iy])
  	}

  	# try to find a smaller value in the immediate neighborhood
  	# make our current position the square with the minimum value.
  	# if a minimum value other than our own current value cannot be
  	# found then we are at a minimum.
  	#
  	# search the neighborhood; three different cases: inner element, corner element, side element
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
		# move to the square with the smaller value, i.e., call compute.centroid on this new square
		# Note: the RETURNED x-y coords in the centroid.x and centroid.y matrix at the current location
		# return the RETURNED x-y coordinates
		if (min.x != ix || min.y != iy)
    {
			r.val <- compute.centroid(min.x,min.y)

			# if explicit is set show the exact connected component
			# otherwise construct a connected componenent where all
			# nodes are connected to a centrol node
			if (explicit)
      {
				centroid.x[ix,iy] <<- min.x
				centroid.y[ix,iy] <<- min.y
				list(bestx=min.x,besty=min.y)
			}
			else
      {
				centroid.x[ix,iy] <<- r.val$bestx
				centroid.y[ix,iy] <<- r.val$besty
				r.val
			}
		}
		#else
		# we have found a minimum
		# note the current x-y in the centroid.x and centroid.y matrix
		# return the current x-y coordinates
		else
    {
			centroid.x[ix,iy] <<- ix
			centroid.y[ix,iy] <<- iy
			list(bestx=ix,besty=iy)
		}
	} # end function compute.centroid
  ###########################################################################

	### iterate over the map and find the centroid for each element
	for (i in 1:xdim)
  {
		for (j in 1:ydim)
    {
			compute.centroid(i,j)
		}
	}

	list(centroid.x=centroid.x,centroid.y=centroid.y)
}

### compute.heat -- compute a heat value map representation of the given distance matrix
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
		stop("compute.heat: heat map can not be computed for a map with a dimension of 1")

  # local function as a shorthand for rowix
	xl <- function(ix,iy)
  {
    #cat("converting (",ix,",",iy,") to row", ix + (iy-1) *xdim,"\n")
		#ix + (iy-1) * x
    rowix(map,ix,iy)
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
  heat <- smooth.2d(as.vector(heat),x=as.matrix(xycoords),nrow=x,ncol=y,surface=FALSE,theta=2)

	heat
}

### df.var.test -- a function that applies the F-test testing the ratio
#                  of the variances of the two data frames
# parameters:
# - df1,df2 - data frames with the same number of columns
# - conf - confidence level for the F-test (default .95)

df.var.test <- function(df1,df2,conf = .95)
{
	if (length(df1) != length(df2))
        stop("df.var.test: cannot compare variances of data frames")

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
	list(ratio=var.ratio.v,conf.int.lo=var.confintlo.v,conf.int.hi=var.confinthi.v)
}

### df.mean.test -- a function that applies the t-test testing the difference
#                   of the means of the two data frames
# parameters:
# - df1,df2 - data frames with the same number of columns
# - conf - confidence level for the t-test (default .95)

df.mean.test <- function(df1,df2,conf = .95)
{
	if (ncol(df1) != ncol(df2))
        stop("df.mean.test: cannot compare means of data frames")

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
	list(diff=mean.diff.v,conf.int.lo=mean.confintlo.v,conf.int.hi=mean.confinthi.v)
}

### vsom.r - vectorized, unoptimized version of the stochastic SOM training algorithm written entirely in R
vsom.r <- function(data,xdim,ydim,alpha,train)
{
    ### some constants
    dr <- nrow(data)
    dc <- ncol(data)
    nr <- xdim*ydim
    nc <- dc # dim of data and neurons is the same

    ### build and initialize the matrix holding the neurons
    cells <- nr * nc        # no. of neurons times number of data dimensions
    v <- runif(cells,-1,1)  # vector with small init values for all neurons
    # NOTE: each row represents a neuron, each column represents a dimension.
    neurons <- matrix(v,nrow=nr,ncol=nc)  # rearrange the vector as matrix

    ### compute the initial neighborhood size and step
    #nsize <- ceiling(sqrt(xdim^2 + ydim^2))
    nsize <- max(xdim,ydim) + 1
    nsize.step <- ceiling(train/nsize)
    step.counter <- 0 # counts the number of epochs per nsize.step

    # convert a 1D rowindex into a 2D map coordinate
    coord2D <- function(rowix)
    {
        x <- (rowix-1) %% xdim + 1
        y <- (rowix-1) %/% xdim + 1

        c(x,y)
    }

    # constants for the Gamma function
    m <- c(1:nr)                     # a vector with all neuron 1D addresses
    m2Ds <- matrix(coord2D(m),nr,2)  # x-y coordinate of ith neuron: m2Ds[i,] = c(xi,yi)

    ### neighborhood function
    Gamma <- function(c)
    {
        c2D <- m2Ds[c,]                          # lookup the 2D map coordinate for c
        c2Ds <- rep(1,nr) %o% c2D                # a matrix with each row equal to c2D
        d <- sqrt((c2Ds - m2Ds)^2 %*% c(1,1))    # distance vector of each neuron from c in terms of map coords!
        hood <- ifelse(d < nsize*1.5,alpha,0.0)  # if m on the grid is in neigh then alpha else 0.0

        as.vector(hood)
    }

    ### training ###
    ### the epochs loop
    for (epoch in 1:train)
    {
        # hood size decreases in disrete nsize.steps
        step.counter <- step.counter + 1
        if (step.counter == nsize.step)
        {
            step.counter <- 0
            nsize <- nsize - 1
        }

        # create a sample training vector
        ix <- sample(1:dr,1)
        xk <- as.numeric(data[ix,])

        ### competitive step
        xk.m <- rep(1,nr) %o% xk
        diff <- neurons - xk.m
        squ <- diff * diff
        s <- squ %*% rep(1,nc)
        o <- order(s)
        c <- o[1]

        ### update step
        gamma.m <- Gamma(c) %o% rep(1,nc)
        neurons <- neurons - diff * gamma.m
    }

    neurons
}

### vsom.f - vectorized and optimized version of the stochastic SOM training algorithm written in Fortran90
vsom.f <- function(data,xdim,ydim,alpha,train)
{
    ### some constants
    dr <- nrow(data)
    dc <- ncol(data)
    nr <- xdim*ydim
    nc <- dc # dim of data and neurons is the same

    ### build and initialize the matrix holding the neurons
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
                       PACKAGE="popsom")

    # unpack the structure and list in result[1]
    v <- result[1]
    neurons <- matrix(v[[1]],nrow=nr,ncol=nc,byrow=FALSE)  # rearrange the result vector as matrix

    neurons
}

# batchsom.private - a wrapper around the batch training algorithm from 'class'
batchsom.private <- function(data,grid,min.radius,max.radius,train,init,radius.type)
{
    set.seed(10)

    initt <- pmatch(init, c("random","sample"))
    radius.type <- pmatch(radius.type,c("gaussian","bubble"))

    data <- as.matrix(data)
    nd <- nrow(data)
    ng <- nrow(grid$pts)

    xdim<- grid$xdim
    ydim<- grid$ydim

    maxit <- ceiling(train/nd)

    if(initt == 1){
        init <- matrix(NA, grid$xdim * grid$ydim, dim(data)[2])
        mi <- apply(data, 2, min)
        ma <- apply(data, 2, max)
        for (i in 1:(xdim*ydim)){
            init[i,] <- mi + (ma - mi) * runif(ncol(init))
        }
    }
    else if (initt == 2){
        init <- data[sample(1L:nd, ng, replace = FALSE), , drop = FALSE]
    }

    nhbrdist <- as.matrix(dist(grid$pts))
    radii<- seq(max.radius,min.radius,len=maxit)

    for(i in 1:maxit)
    {
        cl <- as.numeric(knn1(init, data, 1L:ng))
        if(radius.type == 1)
            A <- exp(-nhbrdist/(2*radii[i]))[,cl]
        else if (radius.type == 2)
            A <- (nhbrdist <= radii[i])[,cl]
        ind <- rowSums(A) > 0
        init[ind, ] <- A[ind, ] %*% data / rowSums(A)[ind]
    }

    list(classif=cl,codes=init,grid=grid)
}


### get.unique.centroids -- a function that computes a list of unique centroids from
#                           a matrix of centroid locations.
#
# parameters:
# - map is an object of type 'map'
# - centroids - a matrix of the centroid locations in the map
get.unique.centroids <- function(map, centroids){
  # get the dimensions of the map
  xdim <- map$xdim
  ydim <- map$ydim
  xlist <- c()
  ylist <- c()
  x.centroid <- centroids$centroid.x
  y.centroid <- centroids$centroid.y

  for(ix in 1:xdim){
    for(iy in 1:ydim) {
      cx <- x.centroid[ix, iy]
      cy <- y.centroid[ix, iy]

      # Check if the x or y of the current centroid is not in the list and if not
      # append both the x and y coordinates to the respective lists
      if(!(cx %in% xlist) || !(cy %in% ylist)){
        xlist <- c(xlist, cx)
        ylist <- c(ylist, cy)
      }
    }
  }

  # return a list of unique centroid positions
  list(position.x=xlist, position.y=ylist)
}


### list.clusters -- A function to get the clusters as a list of lists.
#
# parameters:
# - map is an object of type 'map'
# - centroids - a matrix of the centroid locations in the map
# - unique.centroids - a list of unique centroid locations
# - umat - a unified distance matrix
list.clusters <- function(map,centroids,unique.centroids,umat){
  centroids.x.positions <- unique.centroids$position.x
  centroids.y.positions <- unique.centroids$position.y
  componentx <- centroids$centroid.x
  componenty <- centroids$centroid.y
  cluster_list <- list()

  for(i in 1:length(centroids.x.positions)){
    cx <- centroids.x.positions[i]
    cy <- centroids.y.positions[i]

    # get the clusters associated with a unique centroid and store it in a list
    cluster_list[i] <- list.from.centroid(map, cx, cy, centroids, umat)
  }
  cluster_list
}

### list.from.centroid -- A funtion to get all cluster elements
#                         associated to one centroid.
#
# parameters:
# - map is an object of type 'map'
# - x - the x position of a centroid
# - y - the y position of a centroid
# - centroids - a matrix of the centroid locations in the map
# - umat - a unified distance matrix
list.from.centroid <- function(map,x,y,centroids,umat){
  centroid.x <- x
  centroid.y <- y
  sum <- 0
  xdim <- map$xdim
  ydim <- map$ydim
  centroid_weight <- umat[centroid.x, centroid.y]
  cluster_list <- c()
  for(xi in 1:xdim){
    for(yi in 1:ydim){
      cx <- centroids$centroid.x[xi, yi]
      cy <- centroids$centroid.y[xi, yi]

      if(cx == centroid.x && cy == centroid.y){
        cweight <- umat[xi, yi]
        cluster_list <- c(cluster_list, cweight)
      }
    }
  }
  list(cluster_list)
}

avg.homogeneity <- function(map)
{
  if (is.null(map$labels))
  {
    stop("avg.homogeneity: you need to attach labels to the map")
  }

  ### need to make sure the map doesn't have a dimension of 1
  if (map$xdim <= 1 || map$ydim <= 1)
  {
    stop("avg.homogeneity: map dimensions too small")
  }

  x <- map$xdim
  y <- map$ydim
  nobs <- nrow(map$data)
  centroid.labels <- array(data=list(),dim=c(x,y))

 # find the centroid for each neuron on the map
 centroids <- compute.centroids(map)

 ### attach labels to centroids
 # count the labels in each map cell
 for(i in 1:nobs)
 {
   lab <- as.character(map$labels[i,1])
   nix <- map$visual[i]
   c <- coordinate(map,nix)
   ix <- c[1]
   iy <- c[2]
   cx <- centroids$centroid.x[ix,iy]
   cy <- centroids$centroid.y[ix,iy]
   centroid.labels[[cx,cy]] <- append(centroid.labels[[cx,cy]],lab)
 }

 ### compute average homogeneity of the map: h = (1/nobs)*sum_c majority.label_c
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
       #lhh
       print(majority)

       if (nrow(majority) == 1) # only one label
       {
         print(majority[1])
         m.val <- length(label.v)
       }
       else
       {
         print(majority[1,1])
         m.val <- majority[1,2]
       }
       sum.majority <- sum.majority + m.val
     }
   }
 }
 list(homog=sum.majority/nobs, nclust=n.centroids)
}

majority.labels <- function(map)
{
  if (is.null(map$labels))
  {
    stop("majority.labels: you need to attach labels to the map")
  }

  ### need to make sure the map doesn't have a dimension of 1
  if (map$xdim <= 1 || map$ydim <= 1)
  {
    stop("majority.labels: map dimensions too small")
  }

  x <- map$xdim
  y <- map$ydim
  nobs <- nrow(map$data)
  centroid.labels <- array(data=list(),dim=c(x,y))
  majority.labels <- array(data=list(),dim=c(x,y))


  # find the centroid for each neuron on the map
  centroids <- compute.centroids(map)

  ### attach labels to centroids
  # count the labels in each map cell
  for(i in 1:nobs)
  {
   lab <- as.character(map$labels[i,1])
   nix <- map$visual[i]
   c <- coordinate(map,nix)
   ix <- c[1]
   iy <- c[2]
   cx <- centroids$centroid.x[ix,iy]
   cy <- centroids$centroid.y[ix,iy]
   centroid.labels[[cx,cy]] <- append(centroid.labels[[cx,cy]],lab)
  }

  ### attach majority labels to centroids
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
         # just copy the label from the label vector
         majority.labels[[ix,iy]] <- append(majority.labels[[ix,iy]],label.v[1])
       }
       else
       {
         majority.labels[[ix,iy]] <- append(majority.labels[[ix,iy]],levels(majority[1,1])[1])
       }
     }
   }
  }
  majority.labels
}
