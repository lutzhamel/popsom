test_that("map.predict function works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map.build(df,labels,seed=42)
  p<-map.predict(m,df)
  # print(p)

  # spot check the prediction vector for high
  # confidence predictions, one for each class
  expect_equal(p[11,1],"setosa")
  expect_equal(p[56,1],"versicolor")
  expect_equal(p[145,1],"virginica")
})
