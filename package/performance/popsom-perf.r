library(popsom)
library(som)
library(kohonen)
library(MASS)
library(microbenchmark)

# data sets

# iris
data(iris)
df_iris <- subset(iris, select = -Species)

# wines from package kohonen
data(wines)
df_wines <- as.data.frame(scale(wines))

# synthetic data with three clusters
p <- 10
n <- 500
siglarg <- diag(rep(1, p * p), p, p)
means <- c(0, -50, 50)

clusts <- lapply(means, function(mu) mvrnorm(n = n, mu = rep(mu, p), Sigma = siglarg))
df_sim <- as.data.frame(do.call(rbind, clusts))

# benchmark parameters
xdim <- 10
ydim <- 5
r <- sqrt(xdim**2+ydim**2)
iter <- 1000

# benchmark expressions
bm <- function(dat) {
  data.frame(summary(microbenchmark(
    popsom::map.build(dat,
                      xdim=xdim,
   	              ydim=ydim,
	              train = iter,
		      minimal=TRUE),

    som::som(dat,
             xdim = xdim,
	     ydim = ydim,
	     init="random",
             alpha=c(.3,.3),
             alphaType="linear",
             neigh="bubble",
             topol="rect",
             radius=c(r,r),
             rlen=c(1,iter)),

     kohonen::som(as.matrix(dat),
                  rlen=iter,
	          grid=somgrid(xdim=xdim,ydim=ydim,topo="rectangular",neighbourhood.fct="bubble")),

     times=3)))
}

# run the benchmarks

print("times are reported in milliseconds")

print("### Iris ###")
bmr <- bm(df_iris)
d <- data.frame(cbind(c("popsom","som","kohonen"),format(bmr[["mean"]],digits=2)))
names(d)<-c("package","mean time")
print(d)

print("### Wines ###")
bmr <- bm(df_wines)
d <- data.frame(cbind(c("popsom","som","kohonen"),bmr[["mean"]]))
names(d)<-c("package","mean time")
print(d)

print("### Sim ###")
bmr <- bm(df_sim)
d <- data.frame(cbind(c("popsom","som","kohonen"),bmr[["mean"]]))
names(d)<-c("package","mean time")
print(d)
