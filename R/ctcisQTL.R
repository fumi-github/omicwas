#' Cell-Type-Specific QTL analysis
#'
#' Cell-Type-Specific QTL analysis
#'
#' The functionality is almost the same as the ridge test
#' in ctassoc.
#' Association is tested between each row of Y and each row of X.
#' Usually, the former will be a methylation/expression marker,
#' and the latter will be a SNP.
#' To cope with the large number of combinations,
#' the testing is limited to pairs whose position is within
#' the difference specified by \code{max.pos.diff}
#' In addition, a uniform hyperparameter of ridge regression is used
#' for all pairs.
#' We use the \link[ridge]{linearRidge} function of the ridge package.
#'
#' @param X Matrix (or vector) of SNP genotypes; SNPs x samples.
#' @param Xpos Vector of the physical position of X
#' @param W Matrix of proportion of cell types; samples x cell types.
#' @param Y Matrix (or vector) of bulk omics measurements; markers x samples.
#' @param Ypos Vector of the physical position of Y
#' @param C Matrix (or vector) of covariates; samples x covariates.
#' X, Xpos, W, Y, Ypos, C should be numeric.
#' @param max.pos.diff Maximum positional difference to compute cis-QTL.
#' Association is computed between a row of X and a row of Y,
#' only when thet are within this limit.
#' Since the limiting is only by position, the function needs to be run
#' separately for each chromosome.
#' @param nPC A hyperparameter that (indirectly) specifies the penalty
#' for ridge regression.  If \code{NULL}, it is inferred from the data.
#' A unique hyperparameter is applied to all SNP-marker pairs.
#' @param outdir Output directory.
#' @param outfile Output file.
#' @param num.cores Number of CPU cores to use.
#' Full and marginal tests are run in serial, thus num.cores is ignored.
#' @param seed Seed for random number generation.
#' @return The estimate, statistic, p.value are written to the specified file.
#' @seealso ctRUV
#' @examples
#' \donttest{
#' data(GSE79262small)
#' X    = GSE79262small$X
#' Xpos = GSE79262small$Xpos
#' W    = GSE79262small$W
#' Y    = GSE79262small$Y
#' Ypos = GSE79262small$Ypos
#' C    = GSE79262small$C
#' Y    = Y[seq(1, 601, 20), ] # for brevity
#' Ypos = Ypos[seq(1, 601, 20)]
#' ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
#' }
#'
#' @importFrom data.table rbindlist
#' @importFrom dplyr as_tibble
#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom ridge linearRidge
#' @importFrom rlang inform
#' @importFrom stats lm median
#' @importFrom utils getFromNamespace setTxtProgressBar txtProgressBar write.table
#' @export
ctcisQTL = function (X, Xpos, W, Y, Ypos, C = NULL,
                     max.pos.diff = 1e6,
                     nPC = NULL,
                     outdir = tempdir(),
                     outfile = "ctcisQTL.out.txt",
                     num.cores = 1,
                     seed = 123) {

  X = .as.matrix(X, d = "horizontal", nam = "X")
  W = .as.matrix(W, d = "vertical", nam = "W")
  Y = .as.matrix(Y, d = "horizontal", nam = "Y")
  if (!is.null(C)) {
    C = .as.matrix(C, d = "vertical", nam = "C")
  }
  .check_input_cisQTL(X, Xpos, W, Y, Ypos, C)
  X = .rowcenteralize(X)

  cl = makeCluster(num.cores)
  on.exit(stopCluster(cl))

  # First regress out the intercepts,
  # because all variables are regularized in the ridge package.
  if (is.null(C)) {
    tYadjW = lm(y ~ x,
                data = list(y = t(Y), x = W))$residuals
  } else {
    tYadjW = lm(y ~ x,
                data = list(y = t(Y), x = cbind(W, C)))$residuals
  }
  rm(Y)
  gc()
  J = ncol(tYadjW)

  if (is.null(nPC)) {
    inform("Scanning hyperparameters ...")
    result = c()
    set.seed(seed)
    if (J < 50) {
      Jsample = 1:J
    } else {
      Jsample = sample(J, 50)
    }
    imax = length(Jsample)
    pb = txtProgressBar(max = imax, style = 3)
    for (i in 1:imax) {
      j = Jsample[i]
      Kcisindex = (abs(Xpos - Ypos[j]) < max.pos.diff) &
        (matrixStats::rowVars(X) > 0)
      if (! any(Kcisindex)) { next() }
      Xcis = X[Kcisindex, , drop = FALSE]
      # choose SNP with strongest correlation
      k = which.max(abs(as.numeric(
        Xcis %*%
          as.matrix(tYadjW[, j], d = "vertical"))))
      XW = W * as.numeric(Xcis[k, ])
      mod = ridge::linearRidge(y ~ 0 + x,
                               data = list(y = tYadjW[, j],
                                           x = XW))
      result = c(result, mod$chosen.nPCs)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    nPC = median(result)
  }
  inform(paste0("Ridge regression lambda will be chosen based on ",
                nPC,
                " principal component(s)."))

  inform("Ridge regression ...")
  inform(paste0("Writing output to ", file.path(outdir, outfile)))
  write.table(matrix(c("response", "term", "celltype",
                       "estimate", "statistic", "p.value"), nrow = 1),
              file.path(outdir, outfile),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  tYadjWff = ff::ff(
    tYadjW,
    dim = dim(tYadjW),
    dimnames = dimnames(tYadjW))
  rm(tYadjW)
  gc()
  pb = txtProgressBar(max = J, style = 3)
  for (j in 1:J) {
    y = tYadjWff[, j]
    Kcisindex = (abs(Xpos - Ypos[j]) < max.pos.diff) &
      (matrixStats::rowVars(X) > 0)
    if (! any(Kcisindex)) { next() }
    Xcis = X[Kcisindex, , drop = FALSE]
    result =
      parApply(
        cl = cl,
        Xcis,
        1,
        function (x, y, W, nPC) {
          XW = W * x
          mod = ridge::linearRidge(y ~ 0 + x,
                                   data = list(y = y, x = XW),
                                   nPCs = nPC)
          # Use summary.ridgeLinear to get unscaled coefficients,
          # because coefficients in mod$coef or pvals(mod)$coef are scaled.
          fun = getFromNamespace("summary.ridgeLinear", "ridge")
          res = data.frame(fun(mod)$summaries[[1]]$coefficients)
          res$celltype = sub("^x", "", rownames(res))
          rownames(res) = NULL
          res = res[, c("celltype",
                        "Estimate", "t.value..scaled.", "Pr...t..")]
          colnames(res)[2:4] = c("estimate", "statistic", "p.value")
          return(res)
        },
        y, W, nPC)
    result = dplyr::as_tibble(data.table::rbindlist(result, idcol="term"))
    result$statistic = sign(result$estimate) * result$statistic
    result = cbind(data.frame(response = colnames(tYadjWff)[j]),
                   result)
    write.table(result,
                file.path(outdir, outfile),
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                append = TRUE)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  ff::delete(tYadjWff)
  gc()
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
