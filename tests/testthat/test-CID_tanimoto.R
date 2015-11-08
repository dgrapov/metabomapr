#tests
library(metabomapr)
library(dplyr) # this is not loaded properly in the test env

context("Test cid to tanimoto similarity pipeline")


results<-structure(c(1, 0.158469945355191, 0.171122994652406, 0.158469945355191, 
                     1, 0.646788990825688, 0.171122994652406, 0.646788990825688, 1
), .Dim = c(3L, 3L), .Dimnames = list(c("51", "440649", "439161"
), c("51", "440649", "439161")))

input<-c("[]","51", "440649","[]", "439161",  NA )

test_that("testing input", {
  expect_equal(CID_tanimoto(as.numeric(input)),results)
  expect_equal(CID_tanimoto(input),results)
  expect_equal(CID_tanimoto(as.factor(input)),results)
  expect_equal(CID_tanimoto(data.frame(input)),results)
  expect_equal(CID_tanimoto(as.character(input)),results)
})


test_that("testing query length", {
  expect_equal(CID_tanimoto(input,query.length=1),results)
  expect_equal(CID_tanimoto(input,query.length=2),results)
})