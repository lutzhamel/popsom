test_that("map.starburst works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map.build(df,labels,seed=42)

  # nothing to compare to but just run
  # the code and make sure it doesn't crash
  map.starburst(m)

  expect_equal(!is.null(m), TRUE)
})
