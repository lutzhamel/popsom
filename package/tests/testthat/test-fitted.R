test_that("map.fitted function works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map.build(df,labels,seed=42)
  v<-map.fitted(m)

  # note: we cannot test for specific values because
  # different random number generators/OSs will generate
  # different behavior.  Here we make sure that at
  # least the returned structure makes sense.
  expect_equal(length(v) == nrow(df), TRUE)
})
