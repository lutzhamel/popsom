test_that("starburst works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map(df,labels,seed=42)

  # nothing to compare to but just run
  # the code and make sure it doesn't crash
  starburst(m)

  expect_equal(!is.null(m), TRUE)
})
