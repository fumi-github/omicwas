#' Cell-Type-Specific QTL analysis
#'
#' Cell-Type-Specific QTL analysis
#'
#' A function for analyses of QTL, such as eQTL, mQTL, pQTL.
#' The statistical test is almost identical to
#' \code{ctassoc(test =  "nls.identity", regularize = "TRUE")}.
#' Association analysis is performed between each row of Y and each row of X.
#' Usually, the former will be a methylation/expression marker,
#' and the latter will be a SNP.
#' To cope with the large number of combinations,
#' the testing is limited to pairs whose position is within
#' the difference specified by \code{max.pos.diff}; i.e., limited to cis QTL.
#' In detail, this function performs linear ridge regression,
#' whereas nls.identity test of ctassoc actually is nonlinear regression
#' but with f = identity as normalizing transformation).
#' In order to speed up computation, in the first round, the parameters \eqn{\alpha_{h j}} and
#' \eqn{\gamma_{j l}} are fit by linear regression,
#' and in the second round, \eqn{\beta_{h j k}} are fit and tested by
#' linear ridge regression (see documentation for \link[omicwas]{ctassoc}).
#'
#' @param X Matrix (or vector) of SNP genotypes; SNPs x samples.
#' @param Xpos Vector of the physical position of X
#' @param W Matrix of cell type composition; samples x cell types.
#' @param Y Matrix (or vector) of bulk omics measurements; markers x samples.
#' @param Ypos Vector of the physical position of Y
#' @param C Matrix (or vector) of covariates; samples x covariates.
#' X, Xpos, W, Y, Ypos, C should be numeric.
#' @param max.pos.diff Maximum positional difference to compute cis-QTL.
#' Association analysis is performed between a row of X and a row of Y,
#' only when they are within this limit.
#' Since the limiting is only by position, the function needs to be run
#' separately for each chromosome.
#' @param outdir Output directory.
#' @param outfile Output file.
#' @return The estimate, statistic, p.value are written to the specified file.
#' @seealso ctassoc
#' @examples
#' \donttest{
#' data(GSE79262small)
#' X    = GSE79262small$X
#' Xpos = GSE79262small$Xpos
#' W    = GSE79262small$W
#' Y    = GSE79262small$Y
#' Ypos = GSE79262small$Ypos
#' C    = GSE79262small$C
#' X    = X[seq(1, 3601, 100), ] # for brevity
#' Xpos = Xpos[seq(1, 3601, 100)]
#' ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
#' }
#'
#' @importFrom rlang inform
#' @importFrom stats lm var
#' @importFrom utils setTxtProgressBar txtProgressBar write.table
#' @export
ctcisQTL = function (X, Xpos, W, Y, Ypos, C = NULL,
                     max.pos.diff = 1e6,
                     outdir = tempdir(),
                     outfile = "ctcisQTL.out.txt") {

  X = .as.matrix(X, d = "horizontal", nam = "X")
  W = .as.matrix(W, d = "vertical", nam = "W")
  Y = .as.matrix(Y, d = "horizontal", nam = "Y")
  if (!is.null(C)) {
    C = .as.matrix(C, d = "vertical", nam = "C")
  }
  .check_input_cisQTL(X, Xpos, W, Y, Ypos, C)
  X = .rowcentralize(X)

  # First regress out the effects by W and C
  # in order to accelerate computation.
  if (is.null(C)) {
    tYadjW = lm(y ~ x,
                data = list(y = t(Y), x = W))$residuals
  } else {
    tYadjW = lm(y ~ x,
                data = list(y = t(Y), x = cbind(W, C)))$residuals
  }
  rm(Y)
  gc()

  inform(paste0("Writing output to ", file.path(outdir, outfile)))
  write.table(matrix(c("term", "response", "celltype",
                       "estimate", "statistic", "p.value"), nrow = 1),
              file.path(outdir, outfile),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  pb = txtProgressBar(max = nrow(X), style = 3)
  for (i in 1:nrow(X)) {
    x = X[i, ]
    if (var(x) == 0) { next() }
    Kcisindex = (abs(Ypos - Xpos[i]) < max.pos.diff)
    if (! any(Kcisindex)) { next() }
    tYadjWcis = tYadjW[, Kcisindex, drop = FALSE]
    XW = W * x
    result = .lmridgeLW76(XW, tYadjWcis)
    result = cbind(data.frame(term = rownames(X)[i]),
                   result)
    write.table(result,
                file.path(outdir, outfile),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  invisible(1)
}

.check_input_cisQTL = function (X, Xpos, W, Y, Ypos, C) {
  if (ncol(Y) != ncol(X)) {
    abort("Error: ncol(Y) must equal ncol(X)")
  }
  if (ncol(Y) != nrow(W)) {
    abort("Error: ncol(Y) must equal nrow(W)")
  }
  if (!is.null(C) && ncol(Y) != nrow(C)) {
    abort("Error: ncol(Y) must equal nrow(C)")
  }
  if (nrow(X) != length(Xpos)) {
    abort("Error: nrow(X) must equal length(Xpos)")
  }
  if (nrow(Y) != length(Ypos)) {
    abort("Error: nrow(Y) must equal length(Ypos)")
  }
  return(0)
}
