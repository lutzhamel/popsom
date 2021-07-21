test_that("starburst works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map(df,labels,seed=42)

  # nothing to compare to but just run
  # the code and make sure it doesn't crash
  marginal(m,1)

  expect_equal(!is.null(m), TRUE)
})
