\name{popsom-internal}
\alias{bootstrap}
\alias{best.match}
\alias{accuracy}
\alias{coordinate}
\alias{rowix}
\alias{map.graphics.set}\alias{map.graphics.reset}\alias{plots.heat}\alias{compute.centroids}\alias{compute.umat}\alias{compute.heat}\alias{df.var.test}
\alias{df.mean.test}
\alias{vsom.r}
\alias{vsom.f}
\alias{batchsom.private}\title{Internal popsom functions}\description{  Internal popsom Functions.}\usage{
bootstrap(map,conf.int,data.df,k,sample.acc.v)
best.match(map,obs,full=FALSE)
accuracy(map,sample,data.ix)
coordinate(map,rowix)
rowix(map,x,y)
map.graphics.set()
map.graphics.reset(par.vector)
plot.heat(map,heat,explicit=FALSE,comp=TRUE)
compute.centroids(map,heat,explicit=FALSE)
compute.umat(map,smoothing=NULL)
compute.heat(map,d.in,smoothing=NULL)
df.var.test(df1,df2,conf = .95)
df.mean.test(df1,df2,conf = .95)
vsom.r(data,xdim,ydim,alpha,train)
vsom.f(data,xdim,ydim,alpha,train)
batchsom.private(data,grid,min.radius,max.radius,train,init,radius.type)
}\details{  These are not to be called directly by the user.}\keyword{internal}