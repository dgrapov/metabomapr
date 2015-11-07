#tests


context("Test cid to tanimoto similarity pipeline")


results<-structure(c(1, 0.158469945355191, 0.171122994652406, 0.158469945355191, 
                     1, 0.646788990825688, 0.171122994652406, 0.646788990825688, 1
), .Dim = c(3L, 3L), .Dimnames = list(c("51", "440649", "439161"
), c("51", "440649", "439161")))

test_that("full conversion", {
  expect_equal(CID_tanimoto(c("[]","51", "440649","[]", "439161",  NA )),results)
  expect_equal(CID_tanimoto(factor(c("[]","51", "440649","[]", "439161",  NA ))),results)
  expect_equal(CID_tanimoto(data.frame(c("[]","51", "440649","[]", "439161",  NA ))),results)
})