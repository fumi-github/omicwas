# Linear ridge regression
# Regularization parameter lambda is determined by
# an unbiased version of (Lawless & Wang, 1976).
.lmridgeLW76 = function (X, tY) {
  x = svd(X)
  svdd = x$d
  svdd2 = svdd^2
  svdu = x$u
  svdv = x$v
  # coefficient of principal component regression WITHOUT ridge regularization
  betaPCRs = diag(1/svdd) %*% t(svdu) %*% tY
  sigma2s = colSums((tY - svdu %*% (t(svdu) %*% tY))^2) /
    (nrow(X) - ncol(X))

  weightedmean_betaPCRsquared_adjusteds =
    colSums((svdd * betaPCRs)^2 /
              matrix(sigma2s, byrow = TRUE, nrow = ncol(X), ncol = length(sigma2s)) -
              1) /
    sum(svdd2)
  lambdas = 1 / weightedmean_betaPCRsquared_adjusteds
  lambdas[weightedmean_betaPCRsquared_adjusteds < 0] = svdd2[1]

  # coefficient of linear regression WITH ridge regularization
  betaridges =
    svdv %*%
    (1 / (1 +
            matrix(lambdas, byrow = TRUE, nrow = ncol(X), ncol = length(lambdas)) /
            matrix(svdd2, byrow = FALSE, nrow =length(svdd2), ncol = ncol(tY))) *
       betaPCRs)

  varcomponents =
    apply(svdv,
          2,
          function (x) {diag(x %*% t(x))})
  SEs = sqrt(
    matrix(sigma2s, byrow = TRUE, nrow = ncol(X), ncol = length(sigma2s)) *
      (varcomponents %*%
         (matrix(svdd, byrow = FALSE, nrow = length(svdd), ncol = ncol(tY)) /
            (matrix(svdd2, byrow = FALSE, nrow = length(svdd2), ncol = ncol(tY)) +
               matrix(lambdas, byrow = TRUE, nrow = ncol(X), ncol = length(lambdas))))^2) )

  rownames(betaridges) = colnames(X)
  betaridges = as.data.frame(t(betaridges))
  betaridges$response = rownames(betaridges)
  betaridges =
    tidyr::pivot_longer(
      betaridges,
      cols = -c("response"),
      names_to = "celltypeterm",
      values_to = "estimate")

  rownames(SEs) = colnames(X)
  colnames(SEs) = colnames(tY)
  SEs = as.data.frame(t(SEs))
  SEs$response = rownames(SEs)
  SEs = tidyr::pivot_longer(
    SEs,
    cols = -c("response"),
    names_to = "celltypeterm",
    values_to = "SE")

  result = dplyr::left_join(
    betaridges,
    SEs,
    by = c("response", "celltypeterm"))

  result$statistic = result$estimate / result$SE
  result$p.value = pt(- abs(result$statistic), df = nrow(X) - ncol(X)) * 2
  result = dplyr::select(result, -c("SE"))
  return(result)
}
