test_that("map.significance function works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map.build(df,labels,seed=42)
  s<-map.significance(m,graphics=FALSE)

  # spot check the significance vector for most
  # significant feature
  expect_equal(0.4 < s[3] && s[3]< 0.7 ,TRUE)

})
