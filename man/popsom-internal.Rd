\name{popsom-internal}
\alias{map.embed.vm}
\alias{map.embed.ks}
\alias{map.normalize}
\alias{bootstrap}
\alias{best.match}
\alias{accuracy}
\alias{coordinate}
\alias{rowix}
\alias{map.graphics.set}\alias{map.graphics.reset}\alias{plots.heat}\alias{compute.centroids}\alias{compute.umat}\alias{compute.heat}\alias{df.var.test}
\alias{df.mean.test}
\alias{vsom.r}
\alias{vsom.f}
\alias{batchsom.private}
\alias{compute.combined.clusters}
\alias{get.unique.centroids}
\alias{distance.from.centroids}
\alias{cluster.spread}
\alias{distance.between.clusters}
\alias{list.clusters}
\alias{list.from.centroid}
\alias{combine.decision}
\alias{swap.centroids}
\alias{new.centroid}\title{Internal popsom functions}\description{  Internal popsom Functions.}\usage{
map.embed.vm(map,conf.int=.95,verb=FALSE)
map.embed.ks(map,conf.int=.95,verb=FALSE) 
map.normalize(x, byrow = TRUE)
bootstrap(map,conf.int,data.df,k,sample.acc.v)
best.match(map,obs,full=FALSE)
accuracy(map,sample,data.ix)
coordinate(map,rowix)
rowix(map,x,y)
map.graphics.set()
map.graphics.reset(par.vector)
plot.heat(map, heat, explicit=FALSE, comp=TRUE, merge=FALSE, merge.range)
compute.centroids(map,heat,explicit=FALSE)
compute.umat(map,smoothing=NULL)
compute.heat(map,d.in,smoothing=NULL)
df.var.test(df1,df2,conf = .95)
df.mean.test(df1,df2,conf = .95)
vsom.r(data,xdim,ydim,alpha,train)
vsom.f(data,xdim,ydim,alpha,train)
batchsom.private(data,grid,min.radius,max.radius,train,init,radius.type)
compute.combined.clusters(map,heat,explicit,range)
get.unique.centroids(map, centroids)
distance.from.centroids(map, centroids, unique.centroids, heat)
cluster.spread(x, y, umat, centroids, map)
distance.between.clusters(map, centroids, unique.centroids, umat)
list.clusters(map,centroids,unique.centroids,umat)
list.from.centroid(map,x,y,centroids,umat)
combine.decision(within_cluster_dist,distance_between_clusters,range)
swap.centroids(map, x1, y1, x2, y2, centroids)
new.centroid(bmat, centroids, unique.centroids, map)
}\details{  These are not to be called directly by the user.}\keyword{internal}