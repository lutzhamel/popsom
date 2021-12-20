library(popsom)
library(som)
library(kohonen)
library(MASS)
library(microbenchmark)

# data sets

print("### iris ###")
data(iris)
df_iris <- subset(iris, select = -Species)
x<-15
y<-10
iter<-100000
m<-map.build(df_iris,xdim=x,ydim=y,train=iter,seed=42)
map.summary(m)

print("### epil ###")
# epil from MASS package
data(epil)
df_epil <- as.data.frame(som::normalize(subset(epil,select=-trt)))
x<-15
y<-10
iter<-100000
m<-map.build(df_epil,xdim=x,ydim=y,train=iter,seed=42)
map.summary(m)

print("### wines ###")
# wines from package kohonen
data(wines)
df_wines <- as.data.frame(som::normalize(wines))
x<-15
y<-10
iter<-100000
m<-map.build(df_wines,xdim=x,ydim=y,train=iter,seed=42)
map.summary(m)

print("### synth data ###")
# synthetic data with three clusters
p <- 10
#n <- 500
n <- 100
siglarg <- diag(rep(1, p * p), p, p)
means <- c(0, -50, 50)
clusts <- lapply(means, function(mu) mvrnorm(n = n, mu = rep(mu, p), Sigma = siglarg))
df_sim <- as.data.frame(do.call(rbind, clusts))
x<-25
y<-20
iter<-1000000
m<-map.build(df_sim,xdim=x,ydim=y,train=iter,seed=42)
map.summary(m)
#map.starburst(m)


# function to benchmark algorithms
bm <- function(dat,xdim,ydim,iter)
{
  r <- sqrt(xdim**2+ydim**2)
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
                  rlen=as.integer(iter/nrow(dat)),
                  radius=r,
                  alpha=c(.3,.001),
                  mode="online",
	          grid=somgrid(xdim=xdim,ydim=ydim,topo="rectangular",neighbourhood.fct="bubble")),

     times=3)))
}

# run the benchmarks

# compute speedup vector
speedup <- function (bmr) {
  popsom.time <- bmr[1,"mean"]
  c(popsom.time/popsom.time,bmr[2,"mean"]/popsom.time,bmr[3,"mean"]/popsom.time)
}

print("### Iris ###")
bmr <- bm(df_iris,xdim=15,ydim=10,iter=100000)
d <- data.frame(cbind(c("popsom","som","kohonen"),
                      bmr[["mean"]],
                      speedup(bmr)))
names(d)<-c("package","mean time","speedup")
print(d)

print("### Epil ###")
bmr <- bm(df_epil,xdim=15,ydim=10,iter=100000)
d <- data.frame(cbind(c("popsom","som","kohonen"),
                      bmr[["mean"]],
                      speedup(bmr)))
names(d)<-c("package","mean time","speedup")
print(d)

print("### Wines ###")
bmr <- bm(df_wines,xdim=15,ydim=10,iter=100000)
d <- data.frame(cbind(c("popsom","som","kohonen"),
                      bmr[["mean"]],
                      speedup(bmr)))
names(d)<-c("package","mean time","speedup")
print(d)

print("### Sim ###")
x <- 25
y <- 20
kohonen_df <- as.data.frame(df_sim[sample(1:nrow(df_sim),x*y,replace=TRUE),])
bmr <- bm(kohonen_df,xdim=x,ydim=y,iter=1000000)
d <- data.frame(cbind(c("popsom","som","kohonen"),
                      bmr[["mean"]],
                      speedup(bmr)))
names(d)<-c("package","mean time","speedup")
print(d)
