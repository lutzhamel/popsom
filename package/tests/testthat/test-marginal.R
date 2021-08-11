test_that("map.marginal works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map.build(df,labels,seed=42)

  # nothing to compare to but just run
  # the code and make sure it doesn't crash
  map.marginal(m,1)

  expect_equal(!is.null(m), TRUE)
})
