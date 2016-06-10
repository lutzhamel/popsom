\name{popsom-internal}
\alias{bootstrap}
\alias{accuracy}
\alias{coordinate}
\alias{rowix}
\alias{map.graphics.set}\alias{map.graphics.reset}\alias{plots.heat}\alias{compute.internal.nodes}\alias{compute.umat}\alias{compute.plane}\alias{compute.heat}\alias{df.var.test}
\alias{df.mean.test}\title{Internal popsom functions}\description{  Internal popsom Functions.}\usage{
bootstrap(map,conf.int,data.df,k,sample.acc.v) 
accuracy(map,sample,data.ix)
coordinate(map,rowix)
rowix(map,x,y)
map.graphics.set()map.graphics.reset(par.vector)compute.internal.nodes(map, heat, explicit = FALSE)compute.umat(map, smoothing = NULL)compute.plane(map, plane,smoothing = NULL)compute.heat(map, d, smoothing = NULL)df.var.test(df1, df2, conf = .95)
df.mean.test(df1,df2,conf = .95)plots.heat(map, heat, explicit = FALSE, comp = TRUE)}\details{  These are not to be called directly by the user.}\keyword{internal}