context("ctcisQTL")
library(omicwas)

test_that("ctcisQTL works", {
  skip_on_cran()
  data(GSE79262small)
  X    = GSE79262small$X
  Xpos = GSE79262small$Xpos
  W    = GSE79262small$W
  Y    = GSE79262small$Y
  Ypos = GSE79262small$Ypos
  C    = GSE79262small$C
  Y    = Y[seq(1, 601, 20), ] # for brevity
  Ypos = Ypos[seq(1, 601, 20)]
  ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
  new = read.table(file.path(tempdir(), "ctcisQTL.out.txt"),
                   header = TRUE)
  old = read.table("../ctcisQTL.out.txt",
                   header = TRUE)
  expect_equal(
    head(new$statistic, 100),
    old$statistic)
})
