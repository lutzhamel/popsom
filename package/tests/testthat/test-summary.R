test_that("map.summary function works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map.build(df,labels,seed=42)

  # run summary and pull out convergence from the report tables
  s<-map.summary(m,verb=FALSE)
  conv<-s$quality.assessments[["convergence"]]

  expect_equal(conv < 0.8,TRUE)
})
