#' Cell-Type-Specific Association Testing
#'
#' @param X Matrix (or vector) of phenotypes (and covariates); samples x phenotype(s).
#' @param W Matrix of proportion of cell types; samples x cell types.
#' @param Y Matrix of bulk omics measurements; probes x samples.
#' X, W, Y should be numeric.
#' @param test Statistical test to apply; either ridge, full or marginal.
#' @param num.cores Number of CPU cores to use.
#' Full and marginal tests are run in serial, thus num.cores is ignored.
#' @param chunk.size The size of job for a CPU core in one batch.
#' If you have many cores but limited memory, and there is a memory failure,
#' decrease this parameter.
#' @return A list with one element, which is named "coefficients".
#' The element gives the estimate, statistic, p.value in tibble format.
#' @seealso ctRUV
#' @examples
#' data(GSE42861small)
#' X = GSE42861small$X
#' Y = GSE42861small$Y
#' W = GSE42861small$W
#' Y = ctRUV(X, W, Y)
#' result = ctassoc(X, W, Y)
#' result$coefficients
#'
#' @importFrom broom tidy
#' @importFrom data.table rbindlist
#' @importFrom dplyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom magrittr %>%
#' @importFrom matrixStats colSds
#' @importFrom parallel makeCluster
#' @importFrom parallel parApply
#' @importFrom parallel stopCluster
#' @importFrom ridge linearRidge
#' @importFrom rlang abort
#' @importFrom rlang inform
#' @importFrom R.utils withTimeout
#' @importFrom tidyr pivot_longer
#' @export
ctassoc = function (X, W, Y, test = "ridge",
                    # alpha = 0,
                    # lower.limit = NULL,
                    # upper.limit = NULL,
                    num.cores = 1,
                    chunk.size = 100) {
  if (!(test %in% c("ridge", "full", "marginal"))) {
    abort('Error: test must be either "glmnet", "full" or "marginal"')
  }
  .check_input(X, W, Y)
  X = data.frame(t(t(X)-colMeans(X)))
  switch(test, "ridge" = {
    .full_assoc(X, W, Y, test,
                num.cores = num.cores,
                chunk.size = chunk.size)
  # }, "glmnet" = {
  #   .full_assoc(X, W, Y, test,
  #               alpha = alpha,
  #               lower.limit = lower.limit,
  #               upper.limit = upper.limit,
  #               num.cores = num.cores,
  #               chunk.size = chunk.size)
  }, "full" = {
    .full_assoc(X, W, Y, test,
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "marginal" = {
    .marginal_assoc(X, W, Y)
  })
}

#' Remove Unwanted Variations prior to applying ctassoc
#'
#' @param X Matrix (or vector) of phenotypes (and covariates); samples x phenotype(s).
#' @param W Matrix of proportion of cell types; samples x cell types.
#' @param Y Matrix of bulk omics measurements; probes x samples.
#' X, W, Y should be numeric.
#' @return Y adjusted for the unwanted variations.
#' First, linear regression of Y for the full model is computed.
#' The first ten principal components of the residual is defined as
#'  the unwanted variations, and subtracted from Y.
#' @seealso ctassoc
#' @export
ctRUV = function (X, W, Y) {
  .check_input(X, W, Y)
  X = data.frame(t(t(X)-colMeans(X)))
  X1W = do.call(cbind, apply(W, 2, function(W_h) {cbind(X, 1) * W_h}))
  X1W = as.matrix(X1W)
  YadjX1W = t(lm(y ~ 0 + x,
                 data = list(y = t(Y), x = X1W))$residuals)
  s = svd(YadjX1W)
  rm(YadjX1W)
  gc()
  # regard top 10 principal components as unwanted variation
  s$d[11:length(s$d)] = 0
  D = diag(s$d)
  Y = Y - s$u %*% D %*% t(s$v)
  rm(s, D)
  gc()
  return(Y)
}

.check_input = function (X, W, Y) {
  if (!is.matrix(Y)) {
    abort("Error: Y must be a matrix.")
  }
  if (is.null(dim(W))) {
    abort("Error: W must be a matrix or dataframe.")
  } else {
    if (ncol(Y) != nrow(W)) {
      abort("Error: ncol(Y) must equal nrow(W)")
    }
  }
  if (is.null(dim(X))) {
    if (ncol(Y) != length(X)) {
      abort("Error: ncol(Y) must equal length(vector X)")
    }
  } else {
    if (ncol(Y) != nrow(X)) {
      abort("Error: ncol(Y) must equal nrow(matrix X)")
    }
  }
  return(0)
}

.marginal_assoc = function (X, W, Y) {
  inform("Linear regression, for each cell type")
  Y = t(Y)
  # Avoid parApply, because each process stores Y in memory.
  result = apply(
    W,
    2,
    function (W_h, X, W, Y) {
      cat(".")
      # DEPRECATED; can be heavily biased
      # Xweighted = cbind(X, 1) * W_h
      # res = broom::tidy(lm(y ~ x,
      #                      data = list(y = Y, x = as.matrix(Xweighted))))
      Xweighted = cbind(W, X * W_h)
      res = broom::tidy(lm(y ~ 0 + x,
                           data = list(y = Y, x = as.matrix(Xweighted))))
      res$term = sub("^x", "", res$term)
      return(list(coefficients=res)) },
    X, W, Y)
  cat("\n")
  names(result) = colnames(W)
  return(result)
}

.full_assoc = function (X, W, Y, test,
                        alpha,
                        lower.limit,
                        upper.limit,
                        num.cores,
                        chunk.size) {
  # alpha: 0 for Ridge regression; 1 for Lasso; inbetween for elastic net
  # lower.limit, upper.limit: lowest and highest expression level
  # (in any cell type).
  # The values bind the range of intercepts for the level in a cell type.
  # These should be 0 and 1 for methylation data.
  # If NULL, automatically set to min(Y), max(Y).
  # To disable the binding, set explicitly to -Inf and Inf.
  # Probes with different bound (eg. methylation and gene expression)
  # should not be mixed in one dataset.

  cl = makeCluster(num.cores)
  on.exit(stopCluster(cl))

  Xoriginal = X
  Woriginal = W

  # maintain irreversibility when combining colnames by "." to XW, X1W
  colnames(X) = gsub("\\.", "_", colnames(X))
  colnames(W) = gsub("\\.", "_", colnames(W))
  XW = do.call(cbind, apply(W, 2, function(W_h) {X * W_h}))
  XW = as.matrix(XW)
  X1W = do.call(cbind, apply(W, 2, function(W_h) {cbind(X, 1) * W_h}))
  X1W = as.matrix(X1W)

  switch(test, "full" = { # --------------------------------
    inform("Linear regression ...")
    result = broom::tidy(lm(y ~ 0 + x,
                            data = list(y = t(Y), x = X1W)))
    result$term = sub("^x", "", result$term)
    result = dplyr::rename(result, celltypeterm = term)

  }, "ridge" = { # -----------------------------------------
    inform("Ridge regression ...")
    # First regress out the intercepts,
    # because all variables are regularized in the ridge package.
    tYadjXW = lm(y ~ 0 + x,
                data = list(y = t(Y), x = as.matrix(W)))$residuals
    result =
      .serialparApply(
        cl = cl,
        num.cores = num.cores,
        chunk.size = chunk.size,
        tYadjXW,
        2,
        function (y, XW) {
          # ridge occaisionaly takes very long time
          tryCatch(
            expr = {
              R.utils::withTimeout(
                {mod = ridge::linearRidge(y ~ 0 + x,
                                          data = list(y = y, x = XW))
                # Use summary.ridgeLinear to get unscaled coefficients,
                # because coefficients in mod$coef or pvals(mod)$coef are scaled.
                fun = getFromNamespace("summary.ridgeLinear", "ridge")
                res = data.frame(fun(mod)$summaries[[mod$chosen.nPCs]]$coefficients)
                res$celltypeterm = sub("^x", "", rownames(res))
                rownames(res) = NULL
                res = res[, c("Estimate", "t.value..scaled.", "Pr...t..",
                              "celltypeterm")]
                colnames(res)[1:3] = c("estimate", "statistic", "p.value")
                return(res)
                },
                timeout = 300) },
            TimeoutException = function (ex)
            { data.frame(estimate     = rep(NA, ncol(XW)),
                         statistic    = rep(NA, ncol(XW)),
                         p.value      = rep(NA, ncol(XW)),
                         celltypeterm = colnames(XW)) }) },
        XW)
    rm(tYadjXW)
    gc()
    result = dplyr::as_tibble(data.table::rbindlist(result, idcol="response"))
    result$statistic = sign(result$estimate) * result$statistic

  }, "glmnet" = { # -----------------------------------------
  # For constant terms, don't apply regularization penalty,
  # but bind by (methylation) level in each cell type.
  constantindices = (ncol(X)+1)*(1:ncol(W))
  penalty.factor = rep(1, ncol(X1W))
  penalty.factor[constantindices] = 0
  if (is.null(lower.limit)) {
    lower.limit = min(Y, na.rm=TRUE)
  }
  lower.limit = min(0, lower.limit) # non-positive for glmnet
  if (is.null(upper.limit)) {
    upper.limit = max(Y, na.rm=TRUE)
  }
  lower.limits = rep(-Inf, ncol(X1W))
  lower.limits[constantindices] = lower.limit
  upper.limits = rep(Inf, ncol(X1W))
  upper.limits[constantindices] = upper.limit

  inform("Fitting hyperparameter ...")
  samplingsize = 500
  if (nrow(Y) < samplingsize) {
    Ysmall = Y
  } else {
    set.seed(123)
    Ysmall = Y[sample(nrow(Y), samplingsize), ]
  }
  opt_lambda_list_small =
    parApply(
      cl,
      Ysmall,
      1,
      function (y, X1W, alpha, penalty.factor, lower.limits, upper.limits) {
        set.seed(123)
        cv_fit = glmnet::cv.glmnet(
          x = X1W,
          y = y,
          alpha = alpha,
          lambda = exp(seq(15, -15, -0.5)), # Is this versatile??
          intercept = 0,
          penalty.factor = penalty.factor,
          lower.limits = lower.limits,
          upper.limits = upper.limits)
        return(cv_fit$lambda.min) },
      X1W, alpha, penalty.factor, lower.limits, upper.limits)
  if (nrow(Y) < samplingsize) {
    opt_lambda_list = opt_lambda_list_small
  } else {
    opt_lambda_list = rep(quantile(opt_lambda_list_small, 0.25), nrow(Y))
    names(opt_lambda_list) = NULL
  }

  inform("Regularized regression ...")
  result =
    parApply(
      cl,
      cbind(opt_lambda_list, Y),
      1,
      function (ly, X1W, alpha, penalty.factor, lower.limits, upper.limits) {
        l = ly[1]
        y = ly[-1]
        set.seed(123)
        mod = glmnet::glmnet(
          x = X1W,
          y = y,
          alpha = alpha,
          lambda = l,
          intercept = 0,
          penalty.factor = penalty.factor,
          lower.limits = lower.limits,
          upper.limits = upper.limits)
        return(mod$beta[, 1]) },
      X1W, alpha, penalty.factor, lower.limits, upper.limits)
  result = data.frame(t(result))
  result$response = rownames(result)
  result = tidyr::pivot_longer(
    result,
    cols = -response,
    names_to = "celltypeterm",
    values_to = "estimate")

  inform("Computing statistical significance ...")
  YSd = colSds(lm(y ~ 0 + x,
                  data = list(y = t(Y), x = as.matrix(W))
                  )$residuals)
  # NOT USED: raw Sds for constants, residual Sds for other terms
  # constantindiceslong =
  #   rep(0:(nrow(Y)-1), each=length(constantindices)) * ncol(X1W) +
  #     constantindices
  # YSd[constantindiceslong] = colSds(Y)[constantindiceslong]
  # result$YSd = YSd
  result = result %>%
    left_join(data.frame(response = rownames(Y),
                         YSd = YSd,
                         stringsAsFactors = FALSE),
              by = "response") %>%
    left_join(data.frame(celltypeterm = colnames(X1W),
                         X1WSd = colSds(X1W),
                         stringsAsFactors = FALSE),
              by="celltypeterm") %>%
    # Note: abs(statistic/sqrt(nrow(X1W))) >> 1 occurs for term=="1"
    dplyr::mutate(statistic = sqrt(nrow(X1W))*estimate*X1WSd/YSd) %>%
    dplyr::mutate(p.value = pnorm(abs(statistic), lower.tail = FALSE)*2) %>%
    dplyr::select(-c('YSd', 'X1WSd'))
  }) # end switch ------------------------------

  inform("Summarizing result ...")
  result$celltype =
    c(colnames(Woriginal), "1")[
      match(sub("\\..*", "", result$celltypeterm),
            c(colnames(W), "1"))]
  result$term =
    c(colnames(Xoriginal), "1")[
      match(sub(".*\\.", "", result$celltypeterm),
            c(colnames(X), "1"))]
  result = dplyr::select(
    result,
    c(response, celltype, term, estimate, statistic, p.value))
  return(list(coefficients=result))
}

.serialparApply = function (cl, num.cores, chunk.size, X, MARGIN, FUN, ...) {
  gc()
  batchsize = num.cores * chunk.size
  nbatches = ceiling(dim(X)[MARGIN]/batchsize)
  result = list()
  pb = txtProgressBar(max = nbatches, style = 3)
  for (i in 0:(nbatches - 1)) {
    result = c(
      result,
      parApply(
        cl = cl,
        X = if (MARGIN == 1) { # only works for matrix
          X[seq(1 + i * batchsize,
                min((i+1) * batchsize, nrow(X))), ]
        } else {
          X[, seq(1 + i * batchsize,
                  min((i+1) * batchsize, ncol(X)))]
        },
        MARGIN = MARGIN,
        FUN = FUN,
        ...))
    setTxtProgressBar(pb, i + 1)
    gc()
  }
  close(pb)
  return(result)
}
