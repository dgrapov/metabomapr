#tests
library(metabomapr)
library(reshape2) # this is not loaded properly in the test env

context("Test adjacency to edge list conversion")



input<-structure(c(1, 0.158469945355191, 0.171122994652406, 0.158469945355191, 
                   1, 0.646788990825688, 0.171122994652406, 0.646788990825688, 1
), .Dim = c(3L, 3L), .Dimnames = list(c("51", "440649", "439161"
), c("51", "440649", "439161")))

#lazy results
res1<-structure(list(source = c(51L, 51L, 440649L), target = c(440649L, 439161L, 439161L), value = structure(1:3, .Label = c("0.158469945355191",  "0.171122994652406", "0.646788990825688", "na"), class = "factor")), .Names = c("source",                                                                                                                                                                                              "target", "value"), row.names = c(4L, 7L, 8L), class = "data.frame")

res2<-structure(list(source = c(51L, 51L, 440649L, 51L, 440649L, 439161L
), target = c(51L, 440649L, 440649L, 439161L, 439161L, 439161L
), value = structure(c(4L, 1L, 4L, 2L, 3L, 4L), .Label = c("0.158469945355191", 
                                                           "0.171122994652406", "0.646788990825688", "1", "na"), class = "factor")), .Names = c("source", 
                                                                                                                                                "target", "value"), row.names = c(1L, 4L, 5L, 7L, 8L, 9L), class = "data.frame")

res3<-structure(list(source = c(51L, 440649L, 439161L, 51L, 440649L, 
                                439161L, 51L, 440649L, 439161L), target = c(51L, 51L, 51L, 440649L, 
                                                                            440649L, 440649L, 439161L, 439161L, 439161L), value = structure(c(4L, 
                                                                                                                                              1L, 2L, 1L, 4L, 3L, 2L, 3L, 4L), .Label = c("0.158469945355191", 
                                                                                                                                                                                          "0.171122994652406", "0.646788990825688", "1"), class = "factor")), .Names = c("source", 
                                                                                                                                                                                                                                                                         "target", "value"), row.names = c(NA, 9L), class = "data.frame")


test_that("test conversions", {
  expect_equal(adjacency_edgeList(input,symmetric=TRUE,diagonal=FALSE),res1)
  expect_equal(adjacency_edgeList(input,symmetric=TRUE,diagonal=TRUE),res2)
  expect_equal(adjacency_edgeList(input,symmetric=FALSE,diagonal=TRUE),res3)
})


input2<-structure(c(1, 0.158469945355191, NA, 0.158469945355191, 
                          1, 0.646788990825688, 0.171122994652406, NA, 1
), .Dim = c(3L, 3L), .Dimnames = list(c("51", "440649", "439161"
), c("51", "440649", "439161")))

res4<-structure(list(source = c(51L, 51L, 440649L), target = c(440649L, 
                                                               439161L, 439161L), value = structure(c(1L, 2L, NA), .Label = c("0.158469945355191", 
                                                                                                                              "0.171122994652406", "na", "nna"), class = "factor")), .Names = c("source", 
                                                                                                                                                                                                "target", "value"), row.names = c(4L, 7L, 8L), class = "data.frame")

res5<-structure(list(source = c(440649L, 439161L, 51L, 439161L, 51L, 
                                440649L), target = c(51L, 51L, 440649L, 440649L, 439161L, 439161L
                                ), value = structure(c(1L, NA, 1L, 3L, 2L, NA), .Label = c("0.158469945355191", 
                                                                                           "0.171122994652406", "0.646788990825688", "na", "nna"), class = "factor")), .Names = c("source", 
                                                                                                                                                                                  "target", "value"), row.names = c(2L, 3L, 4L, 6L, 7L, 8L), class = "data.frame")

test_that("testing missing value handling", {
  expect_equal(adjacency_edgeList(input2),res4)
  expect_equal(adjacency_edgeList(input2,symmetric=FALSE),res5)
})