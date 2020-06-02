context("omicwas")
library(omicwas)

test_that("ctassoc works for various test options", {
  skip_on_cran()
  data(GSE42861small)
  X = GSE42861small$X
  Y = GSE42861small$Y
  W = GSE42861small$W
  C = GSE42861small$C
  load("../GSE42861smallresult.RData")
  load("../GSE42861smallresult2.RData")
  expect_equal(
    head(ctRUV(X = X, W = W, Y = Y, C = C), 100),
    resRUV)
  expect_equal(
    head(ctassoc(X = X, W = W, Y = Y, C = C,
                 test = "reducedrankridge")$coefficients$statistic, 100),
    resreducedrankridge$coefficients$statistic)
  expect_equal(
    head(ctassoc(X = X, W = W, Y = Y[1:50, ], C = C,
                 test = "nls.identity",
                 regularize = TRUE,
                 chunk.size = 10)$coefficients$statistic, 100),
    resnlsidentity$coefficients$statistic)
  expect_equal(
    head(ctassoc(X = X, W = W, Y = Y[1:50, ], C = C,
                 test = "nls.logit",
                 regularize = TRUE,
                 chunk.size = 10)$coefficients$statistic, 100),
    resnlslogit$coefficients$statistic)
  expect_equal(
    head(ctassoc(X = X, W = W, Y = Y, C = C,
                 test = "full")$coefficients$statistic, 100),
    resfull$coefficients$statistic)
  expect_equal(
    head(ctassoc(X = X, W = W, Y = Y, C = C,
                 test = "marginal")$CD4.$coefficients$statistic, 100),
    resmarginal$CD4.$coefficients$statistic)
})

test_that("GTEx dataset works for various test options", {
  skip_on_cran()
  data(GTExsmall)
  X = GTExsmall$X
  Y = GTExsmall$Y + 1
  W = GTExsmall$W
  C = GTExsmall$C
  load("../GTExsmallresult.RData")
  expect_equal(
    head(ctassoc(X = X, W = W, Y = Y[1:50, ], C = C,
                 test = "nls.log",
                 regularize = TRUE,
                 chunk.size = 10)$coefficients$statistic, 100),
    resnlslog$coefficients$statistic)
})
