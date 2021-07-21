test_that("predict function works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map(df,labels,seed=42)
  s<-significance(m,graphics=FALSE)

  # spot check the significance vector for most
  # significant feature
  expect_equal(0.4 < s[3] && s[3]< 0.7 ,TRUE)

})
