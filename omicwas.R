library(broom)
library(dplyr)
library(glmnet)
library(parallel)
library(rlang)

.check_input = function (X, W, Y) {
  # X: phenotypes (and covariates); samples x phenotype(s)
  # W: proportion of cell types; samples x cell types
  # Y: bulk omics measurements; probes x samples
  # X, W, Y should be numeric
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

omicassoc = function (X, W, Y, test,
                      alpha = 0,
                      lower.limit = NULL,
                      upper.limit = NULL,
                      num_cores = 1) {
  .check_input(X,W,Y)
  if (!(test %in% c("regularize", "full", "marginal"))) {
    abort('Error: test must be either "regularize", "full" or "marginal"')
  }
  switch(test, "regularize" = {
    .full_assoc(X, W, Y, test,
                alpha = alpha,
                lower.limit = lower.limit,
                upper.limit = upper.limit,
                num_cores = num_cores)
  }, "full" = {
    .full_assoc(X, W, Y, test)
  }, "marginal" = {
    .marginal_assoc(X, W, Y,
                    num_cores = num_cores)
  })
}


.marginal_assoc = function (X, W, Y, num_cores) {
  cl = makeCluster(num_cores)
  on.exit(stopCluster(cl))

  X = data.frame(t(t(X)-colMeans(X)))
  
  Y = t(Y)
  result = parApply(cl,
                    W,
                 2,
                 function(W_h, X, W, Y) {
                   print("*")

# DEPRECATED; can be heavily biased
                          #   Xweighted = cbind(X,1) * W_h
                          #   res = broom::tidy(lm(y ~ x,
                          #                   data=list(y=Y, x=as.matrix(Xweighted))))
                            Xweighted = cbind(W, X * W_h)
                            res = broom::tidy(lm(y ~ 0 + x,
                                          data=list(y=Y, x=as.matrix(Xweighted))))
                   res$term = sub("^x", "", res$term)
                   return(res) },
                 X, W, Y)
  names(result) = colnames(W)
  return(result)
}


.full_assoc = function (X, W, Y, test,
                        alpha,
                        lower.limit,
                        upper.limit,
                        num_cores) {
  # alpha: 0 for Ridge regression; 1 for Lasso; inbetween for elastic net
  # lower.limit, upper.limit: lowest and highest expression level
  # (in any cell type).  If NULL, set as min(Y), max(Y).
  # This should be 0 and 1 for methylation data.
  # Probes with different bound (eg. methylation and gene expression)
  # should not be mixed in one dataset.

  cl = makeCluster(num_cores)
  on.exit(stopCluster(cl))

  Xoriginal = X
  Woriginal = W
  X = data.frame(t(t(X)-colMeans(X)))
  
  # Maintain irreversibility when combining to WX by "."
  colnames(X) = gsub("\\.", "_", colnames(X))
  colnames(W) = gsub("\\.", "_", colnames(W))
  WX = do.call(cbind, apply(W, 2, function(W_h) {cbind(X,1) * W_h}))
  WX = as.matrix(WX)

  switch(test, "regression" = { # ------------------------------
    result = broom::tidy(lm(y ~ 0 + x,
                            data = list(y = t(Y), x = WX)))
    result$term = sub("^x", "", result$term)
    result = rename(result, celltypeterm = term)
    
  },"regularize" = { # -----------------------------------------
  # For constant terms, don't apply regularization penalty,
  # but bind by (methylation) level in each cell type.
  constantindices = (ncol(X)+1)*(1:ncol(W))
  penalty.factor = rep(1, ncol(WX))
  penalty.factor[constantindices] = 0
  if (is.null(lower.limit)) {
    lower.limit = min(Y, na.rm=TRUE)
  }
  lower.limit = min(0, lower.limit) # non-positive for glmnet
  if (is.null(upper.limit)) {
    upper.limit = max(Y, na.rm=TRUE)
  }
  lower.limits = rep(-Inf, ncol(WX))
  lower.limits[constantindices] = lower.limit
  upper.limits = rep(Inf, ncol(WX))
  upper.limits[constantindices] = upper.limit
  
  inform("Inferring hyperparameter ...")
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
      function (y, WX, alpha, penalty.factor, lower.limits, upper.limits) {
        set.seed(123)
        cv_fit = glmnet::cv.glmnet(
          x = WX,
          y = y,
          alpha = alpha,
          intercept = 0,
          penalty.factor = penalty.factor,
          lower.limits = lower.limits,
          upper.limits = upper.limits,
          lambda = exp(seq(15, -15, -0.5))) # Is this versatile??
        return(cv_fit$lambda.min) },
      WX, alpha, penalty.factor, lower.limits, upper.limits)
  if (nrow(Y) < samplingsize) {
    opt_lambda_list = opt_lambda_list_small
  } else {
    opt_lambda_list = rep(quantile(opt_lambda_list_small, 0.25), nrow(Y))
    names(opt_lambda_list) = NULL
  }
  
inform("Regularized regression ...")
result =
  parApply(cl,
  cbind(opt_lambda_list, Y),
  1,
  function (ly, WX, alpha, penalty.factor, lower.limits, upper.limits) {
    l = ly[1]
    y = ly[-1]
    set.seed(123)
    ridge.mod = glmnet::glmnet(x=WX,
                       y=y,
                       alpha=alpha,
                       intercept=0,
                       penalty.factor=penalty.factor,
                       lower.limits=lower.limits,
                       upper.limits=upper.limits,
                       lambda=l)
    return(ridge.mod$beta[,1])
  },
  WX, alpha, penalty.factor, lower.limits, upper.limits)
result = data.frame(t(result))
result$response = rownames(result)
result = tidyr::pivot_longer(result,
                    cols=-response,
                    names_to="celltypeterm",
                    values_to="estimate")
#YSd=rowSds(Y)
YSd=colSds(lm(y ~ 0 + x, data=list(y=t(Y), x=as.matrix(W)))$residuals)
result = result %>%
  left_join(data.frame(response=rownames(Y),
                       YSd = YSd,
                       stringsAsFactors=FALSE),
            by="response") %>%
  left_join(data.frame(celltypeterm=colnames(WX),
                       WXSd=colSds(WX),
                       stringsAsFactors=FALSE),
            by="celltypeterm") %>%
  mutate(statistic=sqrt(nrow(WX))*estimate*WXSd/YSd) %>% # abs(statistic/sqrt(nrow(WX))) >> 1 occurs for term=="1"
  mutate(p.value=pnorm(abs(statistic), lower.tail=FALSE)*2) %>%
  dplyr::select(-c('YSd', 'WXSd'))

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
result = dplyr::select(result, c(response, celltype, term, estimate, statistic, p.value))
return(result)
}
