#' Cell-Type-Specific Association Testing
#'
#' Cell-Type-Specific Association Testing
#'
#' Let the indexes be
#' \eqn{h} for cell type, \eqn{i} for sample,
#' \eqn{j} for marker (CpG site or gene),
#' \eqn{k} for each trait that has cell-type-specific effect,
#' and \eqn{l} for each trait that has a uniform effect across cell types.
#' The input data are \eqn{X_{i k}}, \eqn{C_{i l}}, \eqn{W_{i h}} and \eqn{Y_{j i}},
#' where \eqn{C_{i l}} can be omitted.
#' \eqn{X_{i k}} and \eqn{C_{i l}} are the values for two types of traits,
#' showing effects that are cell-type-specific or not, respectively.
#' Thus, calling \eqn{X_{i k}} and \eqn{C_{i l}} as "traits" and "covariates"
#' gives a rough idea, but is not strictly correct.
#' \eqn{W_{i h}} represents the cell type composition and
#' \eqn{Y_{j i}} represents the marker level,
#' such as methylation or gene expression.
#' For each tissue sample, the cell type proportion \eqn{W_{i h}}
#' is the proportion of each cell type in the bulk tissue,
#' which is measured or imputed beforehand.
#' The marker level \eqn{Y_{j i}} in bulk tissue is measured and provided as input.
#'
#' The parameters we estimate are
#' the cell-type-specific trait effect \eqn{\beta_{h j k}},
#' the tissue-uniform trait effect \eqn{\gamma_{j l}},
#' and the basal marker level \eqn{\alpha_{h j}} in each cell type.
#'
#' We first describe the conventional linear regression models.
#' For marker \eqn{j} in sample \eqn{i},
#' the maker level specific to cell type \eqn{h} is
#' \deqn{\alpha_{h j} + \sum_k \beta_{h j k} * X_{i k}.}
#' This is a representative value rather than a mean, because we do not model
#' a probability distribution for cell-type-specific expression.
#' The bulk tissue marker level is the average weighted by \eqn{W_{i h}},
#' \deqn{\mu_{j i} = \sum_h W_{i h} [ \alpha_{h j} + \sum_k \beta_{h j k} * X_{i k} ] +
#'                   \sum_l \gamma_{j l} C_{i l}.}
#' The statistical model is
#' \deqn{Y_{j i} = \mu_{j i} + \epsilon_{j i},}
#' \deqn{\epsilon_{j i} ~ N(0, \sigma^2_j).}
#' The error of the marker level is is noramlly distributed with variance
#' \eqn{\sigma^2_j}, independently among samples.
#'
#' The \code{full} model is the linear regression
#' \deqn{Y_{j i} = (\sum_h \alpha_{h j} * W_{i h}) +
#'                 (\sum_{h k} \beta_{h j k} * W_{i h} * X_{i k}) +
#'                 (\sum_l \gamma_{j l} * C_{i l}) +
#'                 error.}
#' The \code{marginal} model tests the trait association only in one
#' cell type \eqn{h}, under the linear regression,
#' \deqn{Y_{j i} = (\sum_{h'} \alpha_{h' j} * W_{i h'}) +
#'                 (\sum_k \beta_{h j k} * W_{i h} * X_{i k}) +
#'                 (\sum_l \gamma_{j l} * C_{i l}) +
#'                 error.}
#'
#' The nonlinear model simultaneously analyze cell type composition in
#' linear scale and differential expression/methylation in log/logit scale.
#' The normalizing function is the natural logarithm \eqn{f} = log for gene
#' expression, and \eqn{f} = logit for methylation. Conventional linear regression
#' can be formulated by defining \eqn{f} as the identity function. The three models
#' are named \code{nls.log}, \code{nls.logit} and \code{nls.identity}.
#' We denote the inverse function of \eqn{f} by \eqn{g}; \eqn{g} = exp for
#' gene expression, and \eqn{g} = logistic for methylation.
#' The mean normalized marker level of marker \eqn{j} in sample \eqn{i} becomes
#' \deqn{\mu_{j i} = f(\sum_h W_{i h} g( \alpha_{h j} + \sum_k \beta_{h j k} * X_{i k} )) +
#'                   \sum_l \gamma_{j l} C_{i l}.}
#' The statistical model is
#' \deqn{f(Y_{j i}) = \mu_{j i} + \epsilon_{j i},}
#' \deqn{\epsilon_{j i} ~ N(0, \sigma^2_j).}
#' The error of the marker level is is noramlly distributed with variance
#' \eqn{\sigma^2_j}, independently among samples.
#'
#' The ridge regression aims to cope with multicollinearity of
#' the interacting terms \eqn{W_{i h} * X_{i k}}.
#' Ridge regression is fit by minimizing the residual sum of squares (RSS) plus
#' \eqn{\lambda \sum_{h k} \beta_{h j k}^2}, where \eqn{\lambda > 0} is the
#' regularization parameter.
#'
#' @param X Matrix (or vector) of traits; samples x traits.
#' @param W Matrix of cell type composition; samples x cell types.
#' @param Y Matrix (or vector) of bulk omics measurements; markers x samples.
#' @param C Matrix (or vector) of covariates; samples x covariates.
#' X, W, Y, C should be numeric.
#' @param test Statistical test to apply; either \code{"full"}, \code{"marginal"},
#' \code{"nls.identity"}, \code{"nls.log"}, \code{"nls.logit"} or \code{"reducedrankridge"}.
#' @param regularize Whether to apply Tikhonov (ie ridge) regularization
#' to \eqn{\beta_{h j k}}.
#' The regularization parameter is chosen automatically according to
#' an unbiased version of (Lawless & Wang, 1976).
#' Effective for \code{nls.*} tests.
#' @param num.cores Number of CPU cores to use.
#' Full and marginal tests are run in serial, thus num.cores is ignored.
#' @param chunk.size The size of job for a CPU core in one batch.
#' If you have many cores but limited memory, and there is a memory failure,
#' decrease num.cores and/or chunk.size.
#' @param seed Seed for random number generation.
#' @return A list with one element, which is named "coefficients".
#' The element gives the estimate, statistic, p.value in tibble format.
#' In order to transform the estimate for \eqn{\alpha_{h j}} to the original scale,
#' apply \code{plogis} for \code{test = nls.logit} and
#' \code{exp} for \code{test = nls.log}.
#' The estimate for \eqn{\beta_{h j k}} by \code{test = nls.log} is
#' the natural logarithm of fold-change, not the log2.
#' If numerical convergence fails, \code{NA} is returned for that marker.
#' @references
#' Lawless, J. F., & Wang, P. (1976). A simulation study of ridge and other
#' regression estimators.
#' Communications in Statistics - Theory and Methods, 5(4), 307–323.
#' \url{https://doi.org/10.1080/03610927608827353}
#' @seealso ctcisQTL
#' @examples
#' \donttest{
#' data(GSE42861small)
#' X = GSE42861small$X
#' W = GSE42861small$W
#' Y = GSE42861small$Y
#' C = GSE42861small$C
#' result = ctassoc(X, W, Y, C = C)
#' result$coefficients
#' }
#'
#' @importFrom broom tidy
#' @importFrom data.table rbindlist
#' @importFrom dplyr as_tibble left_join mutate rename select
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom magrittr %>%
#' @importFrom matrixStats colSds
#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom rlang .data abort inform
#' @importFrom stats coef lm median nls nls.control pnorm plogis pt qlogis quantile residuals sd
#' @importFrom tidyr pivot_longer
#' @importFrom utils getFromNamespace setTxtProgressBar txtProgressBar
#' @export
ctassoc = function (X, W, Y, C = NULL,
                    test = "full",
                    regularize = FALSE,
                    # alpha = 0,
                    # lower.limit = NULL,
                    # upper.limit = NULL,
                    num.cores = 1,
                    chunk.size = 1000,
                    seed = 123) {
  if (!(test %in% c("reducedrankridge", "full", "marginal", "nls.identity", "nls.log", "nls.logit"))) {
    abort('Error: test must be either "reducedrankridge", "full", "marginal", "nls.identity", "nls.log", "nls.logit"')
  }
  X = .as.matrix(X, d = "vertical", nam = "X")
  W = .as.matrix(W, d = "vertical", nam = "W")
  Y = .as.matrix(Y, d = "horizontal", nam = "Y")
  if (!is.null(C)) {
    C = .as.matrix(C, d = "vertical", nam = "C")
  }
  .check_input(X, W, Y, C)
  X = .colcenter(X)
  if (!is.null(C)) {
    C = .colcenter(C)
  }
  switch(test, "reducedrankridge" = {
    .full_assoc(X, W, Y, C,
                test = test,
                num.cores = num.cores,
                chunk.size = chunk.size,
                seed = seed)
  }, "nls.identity" = {
    .full_assoc(X, W, Y, C,
                test = "nls",
                nls.link = "identity",
                regularize = regularize,
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "nls.log" = {
    .full_assoc(X, W, Y, C,
                test = "nls",
                nls.link = "log",
                regularize = regularize,
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "nls.logit" = {
    .full_assoc(X, W, Y, C,
                test = "nls",
                nls.link = "logit",
                regularize = regularize,
                num.cores = num.cores,
                chunk.size = chunk.size)
  # }, "glmnet" = {
  #   .full_assoc(X, W, Y, C,
  #               test = test,
  #               alpha = alpha,
  #               lower.limit = lower.limit,
  #               upper.limit = upper.limit,
  #               num.cores = num.cores,
  #               chunk.size = chunk.size,
  #               seed = seed)
  }, "full" = {
    .full_assoc(X, W, Y, C,
                test = test,
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "marginal" = {
    .marginal_assoc(X, W, Y, C)
  })
}

#' Remove Unwanted Variations prior to applying ctassoc
#'
#' Remove Unwanted Variations prior to applying ctassoc
#'
#' First, for each marker, the full linear model of the \code{ctassoc}
#' function is fitted, and the residual is computed.
#' For the residuals over all markers, the principal components (PCs)
#' are computed.
#' The top PCs are regarded as the unwanted variations,
#' and subtracted from \code{Y}.
#'
#' @param X Matrix (or vector) of traits; samples x traits.
#' @param W Matrix of proportion of cell types; samples x cell types.
#' @param Y Matrix (or vector) of bulk omics measurements; markers x samples.
#' @param C Matrix (or vector) of covariates; samples x covariates.
#' X, W, Y, C should be numeric.
#' @param method \code{"PCA"} or \code{"SVA"}
#' @param nPC Number of PCs to be regarded as unwanted variation.
#' If \code{NULL}, automatically computed by the Auer-Gervini approach.
#' @return Y adjusted for the unwanted variations.
#' @seealso ctassoc
#' @importFrom stats formula lm model.matrix
#' @export
ctRUV = function (X, W, Y, C = NULL,
                  method = "PCA",
                  nPC = NULL) {
  X = .as.matrix(X, d = "vertical", nam = "X")
  W = .as.matrix(W, d = "vertical", nam = "W")
  Y = .as.matrix(Y, d = "horizontal", nam = "Y")
  if (!is.null(C)) {
    C = .as.matrix(C, d = "vertical", nam = "C")
  }
  .check_input(X, W, Y, C)
  X = .colcenter(X)
  if (!is.null(C)) {
    C = .colcenter(C)
  }
  X1W = as.matrix(do.call(cbind, apply(W, 2, function(W_h) {cbind(as.data.frame(X), 1) * W_h})))
  switch(method, "PCA" = {
    if (is.null(C)) {
      YadjX1W = t(lm(y ~ x,
                     data = list(y = t(Y), x = X1W))$residuals)
    } else {
      YadjX1W = t(lm(y ~ x,
                     data = list(y = t(Y), x = cbind(X1W, C)))$residuals)
    }
    s = svd(YadjX1W)
    rm(YadjX1W)
    gc()
    if (is.null(nPC)) {
      # nPC computed by Auer-Gervini approach of PCDimension package
      # Same as below
      # spca = ClassDiscovery::SamplePCA(YadjX1W)
      # agDimension(AuerGervini(spca))
      # agDimension(AuerGervini(spca@variances, dd = dim(YadjX1W)))
      nPC = PCDimension::agDimension(PCDimension::AuerGervini(
        Lambda = (s$d)^2/nrow(s$v),
        dd = c(nrow(s$u), nrow(s$v))))
    }
    inform(paste0("Top ", nPC, " PC(s) are regarded as unwanted variation."))
    if (nPC < length(s$d)) {
      s$d[(nPC + 1):length(s$d)] = 0
    }
    D = diag(s$d)
    Y = Y - s$u %*% D %*% t(s$v)
    rm(s, D)
    gc()
  }, "SVA" = {
    if (is.null(C)) {
      mod = model.matrix(formula(paste(c("~ 0", colnames(X1W)), collapse = " + ")),
                         data = as.data.frame(X1W))
      mod0 = model.matrix(formula(paste(c("~ 0", paste0(colnames(W), ".1")), collapse = " + ")),
                          data = as.data.frame(X1W))
    } else {
      mod = model.matrix(formula(paste(c("~ 0", colnames(X1W), colnames(C)), collapse = " + ")),
                         data = as.data.frame(cbind(X1W, C)))
      mod0 = model.matrix(formula(paste(c("~ 0", paste0(colnames(W), ".1"), colnames(C)), collapse = " + ")),
                          data = as.data.frame(cbind(X1W, C)))
    }
    sv = sva::sva(Y,
                  mod,
                  mod0,
                  vfilter = min(nrow(Y), 1e4))$sv
    cat("\n")
    Y = t(lm(t(Y) ~ sv)$residuals)
    rm(sv)
    gc()
  })
  return(Y)
}

.as.matrix = function (X, d, nam = NULL) {
  if (is.null(dim(X))) {
    switch(d, "vertical" = {
      X = matrix(X, ncol = 1)
      if (!is.null(nam)) {
        colnames(X) = nam
      }
    }, "horizontal" = {
      X = matrix(X, nrow = 1)
      if (!is.null(nam)) {
        rownames(X) = nam
      }
    })
  } else {
    X = as.matrix(X)
  }
  return(X)
}

.colcenter = function (m) {
  m -	matrix(colMeans(m, na.rm = TRUE),
             nrow = nrow(m),
             ncol = ncol(m),
             byrow = TRUE)
}

.rowcentralize = function (m) {
  m -	rowMeans(m, na.rm = TRUE)
}

.check_input = function (X, W, Y, C) {
  if (ncol(Y) != nrow(X)) {
    abort("Error: ncol(Y) must equal nrow(X)")
  }
  if (ncol(Y) != nrow(W)) {
    abort("Error: ncol(Y) must equal nrow(W)")
  }
  if (!is.null(C) && ncol(Y) != nrow(C)) {
    abort("Error: ncol(Y) must equal nrow(C)")
  }
  if (ncol(Y) < ncol(X) * ncol(W)) {
    abort("Error: too few observations; ncol(Y) must be at least ncol(X) * ncol(W)")
  }
  return(0)
}

.marginal_assoc = function (X, W, Y, C) {
  inform("Linear regression, for each cell type")
  Y = t(Y)
  # Avoid parApply, because each process stores Y in memory.
  result = apply(
    W,
    2,
    function (W_h, X, W, Y, C) {
      cat(".")
      # DEPRECATED; can be heavily biased
      # Xweighted = cbind(X, 1) * W_h
      # res = broom::tidy(lm(y ~ x,
      #                      data = list(y = Y, x = as.matrix(Xweighted))))
      if (is.null(C)) {
        Xweighted = cbind(W, X * W_h)
      } else {
        Xweighted = cbind(W, X * W_h, C)
      }
      res = broom::tidy(lm(y ~ 0 + x,
                           data = list(y = Y, x = Xweighted)))
      res$term = sub("^x", "", res$term)
      return(list(coefficients = res)) },
    X, W, Y, C)
  cat("\n")
  names(result) = colnames(W)
  return(result)
}

.full_assoc = function (X, W, Y, C,
                        test,
                        nls.link,
                        regularize = TRUE,
                        alpha,
                        lower.limit,
                        upper.limit,
                        num.cores,
                        chunk.size,
                        seed) {
  # alpha: 0 for Ridge regression; 1 for Lasso; inbetween for elastic net
  # lower.limit, upper.limit: lowest and highest expression level
  # (in any cell type).
  # The values bind the range of intercepts for the level in a cell type.
  # These should be 0 and 1 for methylation data.
  # If NULL, automatically set to min(Y), max(Y).
  # To disable the binding, set explicitly to -Inf and Inf.
  # Markers with different bound (eg. methylation and gene expression)
  # should not be mixed in one dataset.

  cl = makeCluster(num.cores)
  on.exit(stopCluster(cl))

  Xoriginal = X
  Woriginal = W

  # maintain irreversibility when combining colnames by "." to XW, X1W
  colnames(X) = gsub("\\.", "_", colnames(X))
  colnames(W) = gsub("\\.", "_", colnames(W))
  X1W = as.matrix(do.call(cbind, apply(W, 2, function(W_h) {cbind(as.data.frame(X), 1) * W_h})))
  XW = X1W[, -c((ncol(X)+1)*(1:ncol(W)))]
  oneXotimesW =
    as.matrix(do.call(cbind, apply(cbind(1, as.data.frame(X)),
                                   2,
                                   function(X_k) {as.data.frame(W) * X_k})))
  colnames(oneXotimesW) =
    sub('([^.]*)\\.([^.]*)', '\\2.\\1', colnames(oneXotimesW), perl = TRUE)

  switch(test, "full" = { # --------------------------------
    inform("Linear regression ...")
    if (is.null(C)) {
      result = broom::tidy(lm(y ~ 0 + x,
                              data = list(y = t(Y), x = X1W)))
    } else {
      result = broom::tidy(lm(y ~ 0 + x,
                              data = list(y = t(Y), x = cbind(X1W, C))))
    }
    result$term = sub("^x", "", result$term)
    result = dplyr::rename(result, celltypeterm = .data$term)

  }, "reducedrankridge" = { # -----------------------------------------
    inform("Reduced-rank ridge regression ...")
    if (is.null(C)) {
      tYadjW = lm(y ~ x,
                  data = list(y = t(Y), x = W))$residuals
    } else {
      tYadjW = lm(y ~ x,
                  data = list(y = t(Y), x = cbind(W, C)))$residuals
    }
    tYadjW = .colcenter(tYadjW)
    rm(Y)
    gc()
    tYadjW_colSds = matrixStats::colSds(tYadjW)
    tYadjWsc = t(t(tYadjW) / tYadjW_colSds)
    rm(tYadjW)
    gc()
    XW = .colcenter(XW)
    XW_colSds = matrixStats::colSds(XW)
    XWsc = t(t(XW) / XW_colSds)

    if (ncol(tYadjWsc) > 1e4) {
      set.seed(seed)
      tYsmall = tYadjWsc[, sample(ncol(tYadjWsc), 1e4)]
    } else {
      tYsmall = tYadjWsc
    }

    inform("Scanning hyperparameters ...")
    imax = 10
    pb = txtProgressBar(max = imax, style = 3)
    SURE = data.frame()
    for (i in 1:imax) {
      lambda = 10^(i - 8)
      rfit = rrs.fit(
        tYsmall,
        XWsc,
        lambda = lambda)
      result = .rrs.SURE(
        rho_ev = rfit$Ad,
        lambda_ev = (svd(XWsc)$d)^2,
        sigma2 = mean((tYsmall - rfit$fitted)^2),
        My = ncol(tYsmall),
        lambda = lambda)
      SURE = rbind(SURE, result)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    rm(tYsmall)
    gc()

    # ggplot(data = SURE) +
    #   geom_point(aes(x = rank,
    #                  y = lambda,
    #                  color = log(Rhat - min(SURE$Rhat) + 1))) +
    #   scale_y_log10()

    SURE = SURE[which.min(SURE$Rhat), ]
    r = SURE$rank[1]
    l = SURE$lambda[1]
    inform(paste0("reduced rank = ", r))
    inform(paste0("ridge penalty lambda = ", l))

    rfit = rrs.fit(
      tYadjWsc,
      XWsc,
      nrank = r,
      lambda = l)
    estimate = rfit$coef
    rm(rfit)

    inform("Jackknife estimation of standard error ...")
    imax = 20
    groupassign = rep(1:imax, ceiling(nrow(tYadjWsc)/imax))[1:nrow(tYadjWsc)]
    set.seed(seed)
    groupassign = groupassign[sample(length(groupassign))]
    coef_jk_ff = ff::ff(
      0,
      dim = c(imax, dim(estimate)))
    pb = txtProgressBar(max = imax, style = 3)
    for (i in 1:imax) {
      coef_jk_ff[i, , ] = rrs.fit(
        tYadjWsc[groupassign != i, ],
        XWsc[groupassign != i, ],
        nrank = r,
        lambda = l)$coef
      gc()
      setTxtProgressBar(pb, i)
    }
    close(pb)
    coef_jk_mean = ff::ffapply(
      X = coef_jk_ff,
      MARGIN = c(2,3),
      AFUN = mean,
      RETURN = TRUE,
      FF_RETURN = FALSE)
    for (i in 1:imax) {
      coef_jk_ff[i, , ] = (coef_jk_ff[i, , ] - coef_jk_mean)^2 }
    SE = sqrt(((imax - 1)/imax) * ff::ffapply(
      X = coef_jk_ff,
      MARGIN = c(2,3),
      AFUN = sum,
      RETURN = TRUE,
      FF_RETURN = FALSE))
    rownames(SE) = rownames(estimate)
    ff::delete(coef_jk_ff)
    rm(coef_jk_mean)
    gc()

    # scale back
    estimate = estimate / XW_colSds
    estimate = t(t(estimate) * tYadjW_colSds)
    SE = SE / XW_colSds
    SE = t(t(SE) * tYadjW_colSds)

    estimate = as.data.frame(t(estimate))
    estimate$response = colnames(tYadjWsc)
    estimate = tidyr::pivot_longer(
      estimate,
      cols = - .data$response,
      names_to = "celltypeterm",
      values_to = "estimate")
    SE = as.data.frame(t(SE))
    SE$response = colnames(tYadjWsc)
    SE = tidyr::pivot_longer(
      SE,
      cols = - .data$response,
      names_to = "celltypeterm",
      values_to = "SE")
    result = estimate %>%
      left_join(SE, by = c("response", "celltypeterm")) %>%
      dplyr::mutate(statistic = estimate/SE) %>%
      dplyr::mutate(p.value = pnorm(abs(.data$statistic), lower.tail = FALSE)*2) %>%
      dplyr::select(-SE)

    ### NOT USED cca or rgcca

    # mineffect = 5 # among the markers
    # # option remove noise
    # if (min(dim(tYadjW)) > 100) {
    #   s = svd(tYadjW)
    #   # remove small PCs as noise
    #   s$d[s$d < sqrt(sum(s$d^2)*mineffect/ncol(tYadjW))] = 0
    #   D = diag(s$d)
    #   tYadjW = s$u %*% D %*% t(s$v)
    # }

    # #   cca = CCA::rcc(tYadjW, XW, 0.001, 0.001) # memory error macbook 658 x 45172
    # #   cca = mixOmics::rcc(tYadjW, XW, lambda1 = 0.001, lambda2 = 0.001)
    #
    # # TODO determine number of partition and then adjust that the sizes are equal
    # tYadjWlist = list()
    # for (i in 1:5) {
    #   tYadjWlist = c(tYadjWlist,
    #                  list(tYadjW[, seq(floor(ncol(tYadjW)*(i-1)/5)+1,
    #                                    floor(ncol(tYadjW)*i/5) )]))
    # }
    # cca = RGCCA::rgcca(c(list(XW), tYadjWlist),
    #                    ncomp = rep(10, length(tYadjWlist)+1))
    #
    # loading = colSums((cca$scores$xscores)^2) / colSums((cca$xcoef)^2)
    # loading =
    #   sqrt(colSums((cca$scores$xscores)^2) /
    #          (sum(tYadjW^2) / ncol(tYadjW))) /
    #   sqrt(colSums((cca$xcoef)^2))
    #
    # topn = sum(CCP::p.asym(
    #   cca$cor,
    #   nrow(tYadjW),
    #   ncol(tYadjW),
    #   ncol(XW))$p.value < 0.05)
    # inform(paste0(topn, " canonical correlations were significant."))

    # estimate = 0
    # SE = 0
    # for (i in 1:topn) {
    #   estimate = estimate +
    #     svdY$d[i] *
    #     matrix(svdY$u[, i], ncol = 1) %*%
    #     matrix(svdYv_coef[[i]]$estimate, nrow = 1)
    #   SE = SE +
    #     svdY$d[i] *
    #     abs(matrix(svdY$u[, i], ncol = 1)) %*%
    #     matrix(svdYv_coef[[i]]$SE, nrow = 1)
    # }
    # estimate = as.data.frame(estimate)
    # SE = as.data.frame(SE)
    # colnames(estimate) = colnames(SE) = colnames(X1W)

  }, "nls" = { # -----------------------------------------
    inform(paste0("nls.", nls.link, " ..."))
    batchsize = num.cores * chunk.size
    totalsize = nrow(Y)
    nbatches = ceiling(totalsize / batchsize)

    switch(
      nls.link,
      logit = {
        if (min(Y, na.rm = TRUE) < 0 | max(Y, na.rm = TRUE) > 1) {
          abort("Error: for test = nls.logit, values of Y must be between 0 and 1")
        }
        if (min(Y, na.rm = TRUE) == 0 | max(Y, na.rm = TRUE) == 1) {
          Y = 0.998 * Y + 0.001
        }
        Y = qlogis(Y)
      }, log = {
        if (min(Y, na.rm = TRUE) <= 0) {
          abort("Error: for test = nls.log, values of Y must be positive")
        }
        Y = log(Y)
      })
    Yff = ff::ff(
      Y,
      dim = dim(Y),
      dimnames = dimnames(Y))
    rm(Y)
    gc()

    mu = switch(nls.link, "identity" = { # --------------------
      if (is.null(C)) {
        function (X, W, oneXotimesW, alpha, beta, sqrtlambda,
                  gradientalpha = FALSE,
                  gradientwithoutalpha = FALSE) {
          beta = matrix(beta, nrow = ncol(W))
          res =
            c(rowSums(W * (rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))),
              sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(oneXotimesW,
                  cbind(matrix(0, nrow = length(beta), ncol = length(alpha)),
                        diag(rep(sqrtlambda, length(beta)))))
          attr(res, "hessian") = function () {
            rep(
              list(diag(rep(0, ncol(oneXotimesW)))),
              nrow(X))
          }
          if (gradientalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, 1:length(alpha)]
          }
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      } else {
        function (X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda,
                  gradientalpha = FALSE,
                  gradientwithoutalpha = FALSE) {
          beta = matrix(beta, nrow = ncol(W))
          res =
            c(rowSums(W * (rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))) +
                C %*% gamma,
              sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(cbind(oneXotimesW, C),
                  cbind(matrix(0, nrow = length(beta), ncol = length(alpha)),
                        diag(rep(sqrtlambda, length(beta))),
                        matrix(0, nrow = length(beta), ncol = length(gamma))))
          attr(res, "hessian") = function () {
            rep(
              list(diag(rep(0, ncol(oneXotimesW) + ncol(C)))),
              nrow(X))
          }
          if (gradientalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, 1:length(alpha)]
          }
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      }
    }, "log" = { # ----------------------------------------
      if (is.null(C)) {
        function (X, W, oneXotimesW, alpha, beta, sqrtlambda,
                  gradientalpha = FALSE,
                  gradientwithoutalpha = FALSE) {
          beta = matrix(beta, nrow = ncol(W))
          g_i_h = exp(rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))
          g_i_h[is.infinite(g_i_h)] = .Machine$double.xmax
          Wlog =
            W * g_i_h /
            (rowSums(W * g_i_h) %*% t(rep(1, ncol(W))))
          res = c(log(rowSums(W * g_i_h)),
                  sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(
              as.matrix(
                do.call(cbind,
                        apply(cbind(1, X),
                              2,
                              function(X_k) {as.data.frame(Wlog) * X_k}))),
              cbind(
                matrix(0, nrow = length(beta), ncol = length(alpha)),
                diag(rep(sqrtlambda, length(beta)))))
          attr(res, "hessian") = function () {
            mapply(
              function (x, w) {
                list(
                  matrix(x, nrow = ncol(X) + 1) %x%
                    matrix(w, nrow = ncol(Wlog))) },
              as.data.frame(
                apply(
                  cbind(1, X),
                  1,
                  function (x) { x %*% t(x) })),
              as.data.frame(
                apply(
                  Wlog,
                  1,
                  function (x) { diag(x) - x %*% t(x) })))
          }
          if (gradientalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, 1:length(alpha)]
          }
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      } else {
        function (X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda,
                  gradientalpha = FALSE,
                  gradientwithoutalpha = FALSE) {
          beta = matrix(beta, nrow = ncol(W))
          g_i_h = exp(rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))
          g_i_h[is.infinite(g_i_h)] = .Machine$double.xmax
          Wlog =
            W * g_i_h /
            (rowSums(W * g_i_h) %*% t(rep(1, ncol(W))))
          res = c(log(rowSums(W * g_i_h)) + C %*% gamma,
                  sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(
              cbind(
                as.matrix(
                  do.call(cbind,
                          apply(cbind(1, X),
                                2,
                                function(X_k) {as.data.frame(Wlog) * X_k}))),
                C),
              cbind(
                matrix(0, nrow = length(beta), ncol = length(alpha)),
                diag(rep(sqrtlambda, length(beta))),
                matrix(0, nrow = length(beta), ncol = length(gamma))))
          attr(res, "hessian") = function () {
            mapply(
              function (x, w) {
                m =
                  matrix(x, nrow = ncol(X) + 1) %x%
                    matrix(w, nrow = ncol(Wlog))
                m = rbind(m,
                          matrix(0, nrow = length(gamma), ncol = ncol(m)))
                m = cbind(m,
                          matrix(0, nrow = nrow(m), ncol = length(gamma)))
                list(m) },
              as.data.frame(
                apply(
                  cbind(1, X),
                  1,
                  function (x) { x %*% t(x) })),
              as.data.frame(
                apply(
                  Wlog,
                  1,
                  function (x) { diag(x) - x %*% t(x) })))
          }
          if (gradientalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, 1:length(alpha)]
          }
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      }
    }, "logit" = { # ----------------------------------------
      if (is.null(C)) {
        function (X, W, oneXotimesW, alpha, beta, sqrtlambda,
                  gradientalpha = FALSE,
                  gradientwithoutalpha = FALSE) {
          beta = matrix(beta, nrow = ncol(W))
          g_i_h = plogis(rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))
          Wlogit =
            W * g_i_h * (1 - g_i_h) /
            ((rowSums(W * g_i_h) * (1 - rowSums(W * g_i_h))) %*% t(rep(1, ncol(W))))
          res = c(qlogis(rowSums(W * g_i_h)),
                  sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(
              as.matrix(
                do.call(cbind,
                        apply(cbind(1, X),
                              2,
                              function(X_k) {as.data.frame(Wlogit) * X_k}))),
              cbind(
                matrix(0, nrow = length(beta), ncol = length(alpha)),
                diag(rep(sqrtlambda, length(beta)))))
          attr(res, "hessian") = function () {
            mapply(
              function (x, w) {
                list(
                  matrix(x, nrow = ncol(X) + 1) %x%
                    matrix(w, nrow = ncol(Wlogit))) },
              as.data.frame(
                apply(
                  cbind(1, X),
                  1,
                  function (x) { x %*% t(x) })),
              as.data.frame(
                apply(
                  (1 - 2 * g_i_h) * Wlogit,
                  1,
                  diag) -
                  sapply(
                    1 - 2 * rowSums(W * g_i_h),
                    function (x) {
                      matrix(x,
                             nrow = ncol(W),
                             ncol = ncol(W))}) *
                  apply(
                    Wlogit,
                    1,
                    function (x) { x %*% t(x) })))
          }
          if (gradientalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, 1:length(alpha)]
          }
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      } else {
        function (X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda,
                  gradientalpha = FALSE,
                  gradientwithoutalpha = FALSE) {
          beta = matrix(beta, nrow = ncol(W))
          g_i_h = plogis(rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))
          Wlogit =
            W * g_i_h * (1 - g_i_h) /
            ((rowSums(W * g_i_h) * (1 - rowSums(W * g_i_h))) %*% t(rep(1, ncol(W))))
          res = c(qlogis(rowSums(W * g_i_h)) + C %*% gamma,
                  sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(
              cbind(
                as.matrix(
                  do.call(cbind,
                          apply(cbind(1, X),
                                2,
                                function(X_k) {as.data.frame(Wlogit) * X_k}))),
                C),
              cbind(
                matrix(0, nrow = length(beta), ncol = length(alpha)),
                diag(rep(sqrtlambda, length(beta))),
                matrix(0, nrow = length(beta), ncol = length(gamma))))
          attr(res, "hessian") = function () {
            mapply(
              function (x, w) {
                m =
                  matrix(x, nrow = ncol(X) + 1) %x%
                  matrix(w, nrow = ncol(Wlogit))
                m = rbind(m,
                          matrix(0, nrow = length(gamma), ncol = ncol(m)))
                m = cbind(m,
                          matrix(0, nrow = nrow(m), ncol = length(gamma)))
                list(m) },
              as.data.frame(
                apply(
                  cbind(1, X),
                  1,
                  function (x) { x %*% t(x) })),
              as.data.frame(
                apply(
                  (1 - 2 * g_i_h) * Wlogit,
                  1,
                  diag) -
                  sapply(
                    1 - 2 * rowSums(W * g_i_h),
                    function (x) {
                      matrix(x,
                             nrow = ncol(W),
                             ncol = ncol(W))}) *
                  apply(
                    Wlogit,
                    1,
                    function (x) { x %*% t(x) })))
          }
          if (gradientalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, 1:length(alpha)]
          }
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      }
    }
    )

    result = list()
    pb = txtProgressBar(max = nbatches, style = 3)
    for (i in 0:(nbatches - 1)) {
      result = c(
        result,
        parApply(
          cl = cl,
          Yff[seq(1 + i * batchsize,
                  min((i+1) * batchsize, totalsize)), ],
          1,
          function (y, X, W, C, oneXotimesW, mu) {
            start_alpha = rep(median(y, na.rm = TRUE), ncol(W))
            lower_alpha = rep(min(y, na.rm = TRUE), ncol(W))
            upper_alpha = rep(max(y, na.rm = TRUE), ncol(W))
            start_beta  = matrix(0, nrow = ncol(W), ncol = ncol(X))
            if (is.null(C)) {
              start_gamma = NULL
            } else {
              start_gamma = rep(0, ncol(C))
            }

            if (is.null(C)) {
              x =
                svd(attr(
                  mu(X, W, oneXotimesW, start_alpha, start_beta, 0),
                  "gradient")[, seq(1, ncol(W) * ncol(X)) + ncol(W)])
            } else {
              x =
                svd(attr(
                  mu(X, W, C, oneXotimesW, start_alpha, start_beta, start_gamma, 0),
                  "gradient")[, seq(1, ncol(W) * ncol(X)) + ncol(W)])
            }
            svdd = x$d
            sqrtlambdalist = c(
              0,
              exp(seq(log(min(svdd)) - 1,
                      log(max(svdd)) + 1,
                      length.out = 20)))

            my_nlsalpha = function (y, X, W, oneXotimesW, C,
                                    start_alpha, start_beta, start_gamma,
                                    lower_alpha,
                                    upper_alpha,
                                    sqrtlambda) {
              if (is.null(C)) {
                mod_alpha = nls(y ~ mu(X, W, oneXotimesW,
                                       alpha, start_beta, sqrtlambda,
                                       gradientalpha = TRUE),
                                data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                            X = as.matrix(X),
                                            W = W,
                                            oneXotimesW = oneXotimesW),
                                start = list(alpha = start_alpha),
                                lower = lower_alpha,
                                upper = upper_alpha,
                                algorithm = "port",
                                control = nls.control(warnOnly = TRUE))
              } else {
                mod_alpha = nls(y ~ mu(X, W, C, oneXotimesW,
                                       alpha, start_beta, start_gamma, sqrtlambda,
                                       gradientalpha = TRUE),
                                data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                            X = as.matrix(X),
                                            W = W,
                                            C = as.matrix(C),
                                            oneXotimesW = oneXotimesW),
                                start = list(alpha = start_alpha),
                                lower = lower_alpha,
                                upper = upper_alpha,
                                algorithm = "port",
                                control = nls.control(warnOnly = TRUE))
              }
              if (! mod_alpha$convInfo$isConv) {
                return("error")
              } else {
                return(list(start_alpha = coef(mod_alpha),
                            mod_alpha   = mod_alpha))
              }
            }
            my_nlswithoutalpha = function (y, X, W, oneXotimesW, C,
                                           start_alpha, start_beta, start_gamma,
                                           lower_alpha,
                                           upper_alpha,
                                           sqrtlambda) {
              if (is.null(C)) {
                mod = nls(y ~ mu(X, W, oneXotimesW,
                                 start_alpha, beta, sqrtlambda,
                                 gradientwithoutalpha = TRUE),
                          data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                      X = as.matrix(X),
                                      W = W,
                                      oneXotimesW = oneXotimesW),
                          start = list(beta = start_beta),
                          algorithm = "port",
                          control = nls.control(warnOnly = TRUE))
                if (mod$convInfo$isConv) {
                  start_beta  = matrix(coef(mod)[seq(1, ncol(W) * ncol(X))],
                                       nrow = ncol(W), ncol = ncol(X))
                } else {
                  return("error")
                }
              } else {
                mod = nls(y ~ mu(X, W, C, oneXotimesW,
                                 start_alpha, beta, gamma, sqrtlambda,
                                 gradientwithoutalpha = TRUE),
                          data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                      X = as.matrix(X),
                                      W = W,
                                      C = as.matrix(C),
                                      oneXotimesW = oneXotimesW),
                          start = list(beta  = start_beta,
                                       gamma = start_gamma),
                          algorithm = "port",
                          control = nls.control(warnOnly = TRUE))
                if (mod$convInfo$isConv) {
                  start_beta  = matrix(coef(mod)[seq(1, ncol(W) * ncol(X))],
                                       nrow = ncol(W), ncol = ncol(X))
                  start_gamma = coef(mod)[seq(1, ncol(C)) + ncol(W) * ncol(X)]
                } else {
                  return("error")
                }
              }
              # compute projection matrix
              if (is.null(C)) {
                x = attr(mu(X, W, oneXotimesW,
                            start_alpha, start_beta,
                            sqrtlambda),
                         "gradient")
              } else {
                x = attr(mu(X, W, C, oneXotimesW,
                            start_alpha, start_beta, start_gamma,
                            sqrtlambda),
                         "gradient")
              }
              e = try({
                P =
                  x[1:length(y), ] %*%
                  solve(t(x) %*% x) %*%
                  t(x[1:length(y), ])
              })
              if (inherits(e, "try-error")) {
                return("error")
              }
              return(list(start_beta  = start_beta,
                          start_gamma = start_gamma,
                          mod         = mod,
                          P           = P))
            }

            nls_result = my_nlsalpha(y, X, W, oneXotimesW, C,
                                     start_alpha, start_beta, start_gamma,
                                     lower_alpha,
                                     upper_alpha,
                                     sqrtlambda = 0)
            # Discard this marker, if nls convergence fails
            if (is.character(nls_result)) {
              res =
                data.frame(estimate     = NA,
                           statistic    = NA,
                           p.value      = NA,
                           celltypeterm = c(colnames(oneXotimesW), colnames(C)))
              return(res)
            }
            start_alpha = nls_result$start_alpha
            mod_alpha   = nls_result$mod_alpha
            RSS = sum((residuals(mod_alpha)[1:length(y)])^2)
            dof_sigma2 = length(y) - length(start_alpha)
            sigma2 = RSS / dof_sigma2

            # sqrtlambda = 0 or the smallest one that nls converges
            for (sqrtlambda in sqrtlambdalist) {
              nls_result = my_nlswithoutalpha(y, X, W, oneXotimesW, C,
                                              start_alpha, start_beta, start_gamma,
                                              lower_alpha,
                                              upper_alpha,
                                              sqrtlambda)
              if (! is.character(nls_result)) { # not "error"
                break()
              }
            }
            # Discard this marker, if nls convergence fails in all of sqrtlambda's.
            if (is.character(nls_result)) {
              res =
                data.frame(estimate     = NA,
                           statistic    = NA,
                           p.value      = NA,
                           celltypeterm = c(colnames(oneXotimesW), colnames(C)))
              return(res)
            }
            start_beta  = nls_result$start_beta
            start_gamma = nls_result$start_gamma
            mod         = nls_result$mod
            P           = nls_result$P

            if (regularize) {

              if (is.null(C)) {
                x =
                  svd(attr(
                    mu(X, W, oneXotimesW, start_alpha, start_beta, 0),
                    "gradient")[, seq(1, ncol(W) * ncol(X)) + ncol(W)])
              } else {
                x =
                  svd(attr(
                    mu(X, W, C, oneXotimesW, start_alpha, start_beta, start_gamma, 0),
                    "gradient")[, seq(1, ncol(W) * ncol(X)) + ncol(W)])
              }
              svdd = x$d
              svdu = x$u
              svdv = x$v

              # # optimal lambda according to [Hoerl 1975]
              # yattributabletobeta =
              #   mu(X, W, C, oneXotimesW, start_alpha, start_beta, start_gamma, sqrtlambda) -
              #   mu(X, W, C, oneXotimesW, start_alpha, start_beta * 0, start_gamma, sqrtlambda)
              # yattributabletobeta = yattributabletobeta[1:length(y)]
              # (t(svdu[1:length(y), ]) %*% yattributabletobeta) / svdd == betaPCR below
              # mean_betasquared_adjusted =
              #   mean(start_beta^2) - mean(sigma2 / svdd^2) # unbiased
              # if (mean_betasquared_adjusted > 0) {
              #   sqrtlambda = sqrt(sigma2 / mean_betasquared_adjusted)
              # } else {
              #   sqrtlambda = svdd[1]
              # }
              betaPCR = t(svdv) %*% start_beta  # alpha in [Cule and Iorio 2013]
              weightedmean_betaPCRsquared_adjusted =
                sum((svdd * betaPCR)^2 / sigma2 - 1) / sum(svdd^2)
              if (weightedmean_betaPCRsquared_adjusted > 0) {
                sqrtlambda = sqrt(1 / weightedmean_betaPCRsquared_adjusted)
              } else {
                sqrtlambda = svdd[1]
              }

              # # optimal lambda according to [Cule and Iorio 2013]
              # # Simplified such that sigma2 is reused.
              # betaPCR = t(svdv) %*% start_beta  # alpha in [Cule and Iorio 2013]
              # dataPCR = data.frame()
              # for (r in 1:length(betaPCR)) {
              #   mean_betaPCRsquared_adjusted =
              #     mean((betaPCR^2 - sigma2 / svdd^2)[1:r]) # unbiased
              #   if (mean_betaPCRsquared_adjusted > 0) {
              #     lambda = sigma2 / mean_betaPCRsquared_adjusted
              #   } else {
              #     lambda = svdd[1]^2
              #   }
              #   dof = sum(1 / (1 + lambda / svdd^2)^2)
              #   dofres = sum(1 - 1 / (1 + svdd^2 / lambda)^2)
              #   twolnlik = sum(
              #     svdd^2 * betaPCR^2 / sigma2 *
              #       (1 - 1 / (1 + svdd^2 / lambda)^2))
              #   MSE =
              #     sum(1 / (1 + svdd^2 / lambda)^2 * betaPCR^2) +
              #     sum(1 / (1 + lambda / svdd^2)^2 * sigma2 / svdd^2)
              #   dataPCR = rbind(dataPCR,
              #                   data.frame(r = r,
              #                              lambda = lambda,
              #                              dof = dof,
              #                              diff = r - dof,
              #                              dofres = dofres,
              #                              AIC = 2 * dof - twolnlik,
              #                              MSE = MSE))
              # }
              # sqrtlambda = sqrt(dataPCR$lambda[which.min(dataPCR$diff)])

              # in case of nls convergence failure, try from similar values
              sqrtlambdalist = c(
                sqrtlambda,
                sqrtlambdalist[order(abs(sqrtlambdalist - sqrtlambda))])

              # avoid inheriting false conversion of start_beta
              start_beta = start_beta * 0

              # GCVdata = data.frame()
              for (sqrtlambda in sqrtlambdalist) {
                nls_result = my_nlswithoutalpha(y, X, W, oneXotimesW, C,
                                                start_alpha, start_beta, start_gamma,
                                                lower_alpha,
                                                upper_alpha,
                                                sqrtlambda)
                if (! is.character(nls_result)) { # not "error"
                  break()
                }
              }

              start_beta  = nls_result$start_beta
              start_gamma = nls_result$start_gamma
              # mod         = nls_result$mod
              # P           = nls_result$P
              # dof = sum(diag(P))
              # RSS = sum((residuals(mod)[1:length(y)])^2)
              # GCV = length(y) * RSS / (length(y) - dof)^2
              # AIC = 2 * dof +
              #   RSS / sigma2 +
              #   length(y) * log(2 * pi * sigma2)
              # BIC = log(length(y)) * dof +
              #   RSS / sigma2 +
              #   length(y) * log(2 * pi * sigma2)
              # GCVdata = rbind(GCVdata,
              #                 data.frame(
              #                   sqrtlambda = sqrtlambda,
              #                   dof = dof,
              #                   GCV = GCV))

            } # regularize

            res = data.frame(estimate = c(start_beta, start_gamma))
            # Wald test
            # To be accurate, non-exact t-type test [Halawa 2000]
            if (is.null(C)) {
              x = attr(mu(X, W, oneXotimesW,
                          start_alpha, start_beta,
                          sqrtlambda),
                       "gradient")
              xx = attr(mu(X, W, oneXotimesW,
                           start_alpha, start_beta,
                           sqrtlambda),
                        "hessian")()
              r = y - mu(X, W, oneXotimesW,
                         start_alpha, start_beta,
                         sqrtlambda)[1:length(y)]
            } else {
              x = attr(mu(X, W, C, oneXotimesW,
                          start_alpha, start_beta, start_gamma,
                          sqrtlambda),
                       "gradient")
              xx = attr(mu(X, W, C, oneXotimesW,
                           start_alpha, start_beta, start_gamma,
                           sqrtlambda),
                        "hessian")()
              r = y - mu(X, W, C, oneXotimesW,
                         start_alpha, start_beta, start_gamma,
                         sqrtlambda)[1:length(y)]
            }
            x = x[, -(1:length(start_alpha))]
            xx = lapply(xx,
                        function (x) {
                          x[-(1:length(start_alpha)),
                            -(1:length(start_alpha))] })
            sigma2Hstar =
              t(x[1:length(y), ]) %*%
              x[1:length(y), ]
            # sigma2Hstarlambdainv = solve(t(x) %*% x)
            # SE = sqrt(sigma2 * diag(sigma2Hstarlambdainv %*%
            #                           sigma2Hstar %*%
            #                           sigma2Hstarlambdainv))
            z = matrix(
              rowSums(
                mapply(
                  function (x, r) { x * r },
                  xx,
                  as.list(r))),
              nrow = nrow(xx[[1]]))
            e = try({ sigma2Hlambdainv = solve(t(x) %*% x - z) })
            if (inherits(e, "try-error")) {
              res =
                data.frame(estimate     = NA,
                           statistic    = NA,
                           p.value      = NA,
                           celltypeterm = c(colnames(oneXotimesW), colnames(C)))
              return(res)
            }
            SE = sqrt(sigma2 * diag(sigma2Hlambdainv %*%
                                      sigma2Hstar %*%
                                      sigma2Hlambdainv))
            res$statistic = res$estimate / SE
            res$p.value = pt(- abs(res$statistic), df = dof_sigma2) * 2
            res_alpha = summary(mod_alpha)$coefficients[, -2]
            colnames(res_alpha) = c("estimate", "statistic", "p.value")
            res = rbind(res_alpha, res)
            res$celltypeterm = c(colnames(oneXotimesW), colnames(C))
            return(res)
          },
          X, W, C, oneXotimesW, mu))
      setTxtProgressBar(pb, i + 1)
    }
    close(pb)
    ff::delete(Yff)
    gc()
    result = dplyr::as_tibble(data.table::rbindlist(result, idcol="response"))

  }, "glmnet" = { # -----------------------------------------
    if (!is.null(C)) {
      tYadjC = lm(y ~ 0 + x,
                  data = list(y = t(Y), x = C))$residuals
      Y = t(tYadjC)
    }
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
    set.seed(seed)
    Ysmall = Y[sample(nrow(Y), samplingsize), ]
  }
  opt_lambda_list_small =
    parApply(
      cl,
      Ysmall,
      1,
      function (y, X1W, alpha, penalty.factor, lower.limits, upper.limits) {
        set.seed(seed)
        cv_fit = glmnet::cv.glmnet(
          x = X1W,
          y = y,
          alpha = alpha,
          lambda = exp(seq(15, -15, -0.5)), # Is this versatile??
          intercept = 0,
          penalty.factor = penalty.factor,
          lower.limits = lower.limits,
          upper.limits = upper.limits,
          standardize = FALSE)
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
        set.seed(seed)
        mod = glmnet::glmnet(
          x = X1W,
          y = y,
          alpha = alpha,
          lambda = l,
          intercept = 0,
          penalty.factor = penalty.factor,
          lower.limits = lower.limits,
          upper.limits = upper.limits,
          standardize = FALSE)
        return(mod$beta[, 1]) },
      X1W, alpha, penalty.factor, lower.limits, upper.limits)
  result = as.data.frame(t(result))
  result$response = rownames(result)
  result = tidyr::pivot_longer(
    result,
    cols = - .data$response,
    names_to = "celltypeterm",
    values_to = "estimate")

  inform("Computing statistical significance ...")
  YSd = colSds(lm(y ~ 0 + x,
                  data = list(y = t(Y), x = W)
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
    dplyr::mutate(statistic = sqrt(nrow(X1W))*estimate*.data$X1WSd/YSd) %>%
    dplyr::mutate(p.value = pnorm(abs(.data$statistic), lower.tail = FALSE)*2) %>%
    dplyr::select(-c("YSd", "X1WSd"))
  }) # end switch ------------------------------

  inform("Summarizing result ...")
  result$celltype =
    c(colnames(Woriginal), "1")[
      match(sub("\\..*", "", result$celltypeterm),
            c(colnames(W), "1"))]
  result$term =
    c(colnames(Xoriginal), "1", colnames(C))[
      match(sub(".*\\.", "", result$celltypeterm),
            c(colnames(X), "1", colnames(C)))]
  result = dplyr::select(
    result,
    c("response", "celltype", "term", "estimate", "statistic", "p.value"))
  return(list(coefficients = result))
}
