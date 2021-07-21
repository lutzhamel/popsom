test_that("position function works", {
  data(iris)
  df<-subset(iris,select=-Species)
  labels<-subset(iris,select=Species)
  m<-map(df,labels,seed=42)
  p<-position(m,df)
  # print(p)

  # note: we cannot test for specific values because
  # different random number generators/OSs will generate
  # different behavior.  Here we make sure that at
  # least the returned structure makes sense.

  expect_equal(nrow(p)==nrow(df),TRUE)
  expect_equal(ncol(p), 2)
})
