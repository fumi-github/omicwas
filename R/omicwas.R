#' Cell-Type-Specific Association Testing
#'
#' Cell-Type-Specific Association Testing
#'
#' Let the indexes be
#' \eqn{h} for cell type, \eqn{i} for sample,
#' \eqn{j} for marker (CpG site or gene),
#' \eqn{k} for each trait that has cell-type-specific effect,
#' and \eqn{l} for each trait that has bulk tissue effect.
#' The input data are \eqn{X_{i k}}, \eqn{C_{i l}}, \eqn{W_{i h}} and \eqn{Y_{j i}},
#' where \eqn{C_{i l}} can be omitted.
#' \eqn{X_{i k}} and \eqn{C_{i l}} are the two types of traits,
#' showing effects that are cell-type-specific or not, respectively.
#' Thus, calling \eqn{X_{i k}} and \eqn{C_{i l}} as "traits" and "covariates"
#' gives a rough idea, but is not strictly correct.
#' \eqn{W_{i h}} represents the cell type proportion and
#' \eqn{Y_{j i}} represents the marker level,
#' such as methylation or gene expression.
#' For each tissue sample, the cell type proportion \eqn{W_{i h}}
#' is the proportion of each cell type in the bulk tissue,
#' which is measured or imputed beforehand.
#' The marker level \eqn{Y_{j i}} in bulk tissue is measured and provided as input.
#'
#' The cell-type-specific marker level \eqn{Z_{h j i}} is not observed
#' and is treated as a hidden variable.
#' The parameters we estimate are
#' the effect of cell-type-specific traits \eqn{\beta_{h j k}},
#' the effect of non-specific traits \eqn{\gamma_{j l}},
#' and the cell-type-specific basal marker level \eqn{\mu_{h j}}.
#'
#' We assume normal distribution for the cell-type-specific marker level,
#' \deqn{Z_{h j i} ~ N(\mu_{h j} + \sum_k \beta_{h j k} * X_{i k}, \sigma^2_{h j}).}
#' Since the bulk tissue marker level is the sum weighted by \eqn{W_{i h}},
#' \deqn{Y_{j i} ~ N(\sum_h W_{i h} [\mu_{h j} + \sum_k \beta_{h j k} * X_{i k}] +
#'                   \sum_l \gamma_{j l} C_{i l}, \tau^2_j).}
#' Although formally, the variance comprises of components of cell type level
#' and tissue level, we approximate and unify into \eqn{\tau^2_j}.
#'
#' The \code{full} model is the linear regression
#' \deqn{Y_{j i} ~ (\sum_h \mu_{h j} * W_{i h}) +
#'                 (\sum_{h k} \beta_{h j k} * W_{i h} * X_{i k}) +
#'                 (\sum_l \gamma_{j l} * C_{i l}) +
#'                 error.}
#' The \code{ridge} model aims to cope with multicollinearity of
#' the interacting terms \eqn{W_{i h} * X_{i k}}.
#' It first adjusts for \eqn{\mu_{h j}} and \eqn{\gamma_{j l}}
#' by fitting linear regression and taking the residuals.
#' Afterwards, ridge regression is used to fit \eqn{\beta_{h j k}}.
#' We use the \link[ridge]{linearRidge} function of the ridge package.
#' The \code{marginal} model tests the trait association only in one
#' cell type \eqn{h}, under the linear regression,
#' \deqn{Y_{j i} ~ (\sum_{h'} \mu_{h' j} * W_{i h'}) +
#'                 (\sum_k \beta_{h j k} * W_{i h} * X_{i k}) +
#'                 (\sum_l \gamma_{j l} * C_{i l}) +
#'                 error.}
#'
#' @param X Matrix (or vector) of traits; samples x traits.
#' @param W Matrix of proportion of cell types; samples x cell types.
#' @param Y Matrix (or vector) of bulk omics measurements; markers x samples.
#' @param C Matrix (or vector) of covariates; samples x covariates.
#' X, W, Y, C should be numeric.
#' @param test Statistical test to apply; either \code{"reducedrankridge"},
#' \code{"ridge"}, \code{"full"} or \code{"marginal"}.
#' @param num.cores Number of CPU cores to use.
#' Full and marginal tests are run in serial, thus num.cores is ignored.
#' @param chunk.size The size of job for a CPU core in one batch.
#' If you have many cores but limited memory, and there is a memory failure,
#' decrease num.cores and/or chunk.size.
#' @param seed Seed for random number generation.
#' @return A list with one element, which is named "coefficients".
#' The element gives the estimate, statistic, p.value in tibble format.
#' @seealso ctRUV
#' @examples
#' \donttest{
#' data(GSE42861small)
#' X = GSE42861small$X
#' Y = GSE42861small$Y
#' W = GSE42861small$W
#' C = GSE42861small$C
#' Y = ctRUV(X, W, Y, C = C)
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
#' @importFrom ridge linearRidge
#' @importFrom rlang .data abort inform
#' @importFrom R.utils withTimeout
#' @importFrom stats coef lm nls nls.control pnorm plogis pt qlogis quantile residuals sd
#' @importFrom tidyr pivot_longer
#' @importFrom utils getFromNamespace setTxtProgressBar txtProgressBar
#' @export
ctassoc = function (X, W, Y, C = NULL,
                    test = "ridge",
                    # alpha = 0,
                    # lower.limit = NULL,
                    # upper.limit = NULL,
                    num.cores = 1,
                    chunk.size = 1000,
                    seed = 123) {
  if (!(test %in% c("reducedrankridge", "ridge", "full", "marginal", "nls.identity", "nls.log", "nls.logit"))) {
    abort('Error: test must be either "reducedrankridge", "ridge", "full", "marginal", "nls.identity", "nls.log", "nls.logit"')
  }
  X = .as.matrix(X, d = "vertical", nam = "X")
  W = .as.matrix(W, d = "vertical", nam = "W")
  Y = .as.matrix(Y, d = "horizontal", nam = "Y")
  if (!is.null(C)) {
    C = .as.matrix(C, d = "vertical", nam = "C")
  }
  .check_input(X, W, Y, C)
  X = .colcentralize(X)
  if (!is.null(C)) {
    C = .colcentralize(C)
  }
  switch(test, "reducedrankridge" = {
    .full_assoc(X, W, Y, C,
                test = test,
                num.cores = num.cores,
                chunk.size = chunk.size,
                seed = seed)
  }, "ridge" = {
    .full_assoc(X, W, Y, C,
                test = test,
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "nls.identity" = {
    .full_assoc(X, W, Y, C,
                test = "nls",
                nls.link = "identity",
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "nls.log" = {
    .full_assoc(X, W, Y, C,
                test = "nls",
                nls.link = "log",
                num.cores = num.cores,
                chunk.size = chunk.size)
  }, "nls.logit" = {
    .full_assoc(X, W, Y, C,
                test = "nls",
                nls.link = "logit",
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
  X = .colcentralize(X)
  if (!is.null(C)) {
    C = .colcentralize(C)
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

.colcentralize = function (m) {
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
    tYadjW = .colcentralize(tYadjW)
    rm(Y)
    gc()
    tYadjW_colSds = matrixStats::colSds(tYadjW)
    tYadjWsc = t(t(tYadjW) / tYadjW_colSds)
    rm(tYadjW)
    gc()
    XW = .colcentralize(XW)
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

  }, "ridge" = { # -----------------------------------------
    inform("Ridge regression ...")
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
    batchsize = num.cores * chunk.size
    totalsize = ncol(tYadjW)
    nbatches = ceiling(totalsize/batchsize)
    tYadjWff = ff::ff(
      tYadjW,
      dim = dim(tYadjW),
      dimnames = dimnames(tYadjW))
    rm(tYadjW)
    gc()
    result = list()
    pb = txtProgressBar(max = nbatches, style = 3)
    for (i in 0:(nbatches - 1)) {
      result = c(
        result,
        parApply(
          cl = cl,
          tYadjWff[, seq(1 + i * batchsize,
                         min((i+1) * batchsize, totalsize))],
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
              TimeoutException = function (ex) {
                data.frame(estimate     = rep(NA, ncol(XW)),
                           statistic    = rep(NA, ncol(XW)),
                           p.value      = rep(NA, ncol(XW)),
                           celltypeterm = colnames(XW)) }) },
          XW))
      setTxtProgressBar(pb, i + 1)
    }
    close(pb)
    ff::delete(tYadjWff)
    gc()
    result = dplyr::as_tibble(data.table::rbindlist(result, idcol="response"))
    result$statistic = sign(result$estimate) * result$statistic

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
                  gradientwithoutalpha = FALSE) {
          res =
            c(rowSums(W * (rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))),
              sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(oneXotimesW,
                  cbind(matrix(0, nrow = length(beta), ncol = length(alpha)),
                        diag(rep(sqrtlambda, length(beta)))))
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      } else {
        function (X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda,
                  gradientwithoutalpha = FALSE) {
          res =
            c(rowSums(W * (rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))) +
                C %*% gamma,
              sqrtlambda * beta)
          attr(res, "gradient") =
            rbind(cbind(oneXotimesW, C),
                  cbind(matrix(0, nrow = length(beta), ncol = length(alpha)),
                        diag(rep(sqrtlambda, length(beta))),
                        matrix(0, nrow = length(beta), ncol = length(gamma))))
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      }
    }, "log" = { # ----------------------------------------
      if (is.null(C)) {
        function (X, W, oneXotimesW, alpha, beta, sqrtlambda,
                  gradientwithoutalpha = FALSE) {
          res = 0 # TODO

          attr(res, "gradient") = 0

          return(res)
        }
      } else {
        function (X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda,
                  gradientwithoutalpha = FALSE) {
          g_i_h = exp(rep(1, nrow(X)) %*% t(alpha) + X %*% t(beta))
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
          if (gradientwithoutalpha) {
            attr(res, "gradient") = attr(res, "gradient")[, -(1:length(alpha))]
          }
          return(res)
        }
      }
    }, "logit" = { # ----------------------------------------
      if (is.null(C)) {
        function (X, W, oneXotimesW, alpha, beta, sqrtlambda,
                  gradientwithoutalpha = FALSE) {
          res = 0 # TODO

          attr(res, "gradient") = 0

          return(res)
        }
      } else {
        function (X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda,
                  gradientwithoutalpha = FALSE) {
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
            GCVstats = data.frame()
            start_alphas = list()
            start_betas  = list()
            start_gammas = list()
            sigma2 = NA

            start_alpha_0 =
              start_alpha = rep(median(y, na.rm = TRUE), ncol(W))
            lower_alpha   = rep(min(y, na.rm = TRUE), ncol(W))
            upper_alpha   = rep(max(y, na.rm = TRUE), ncol(W))
            start_beta_0 =
              start_beta  = matrix(   0, nrow = ncol(W), ncol = ncol(X))
            lower_beta    = matrix(-Inf, nrow = ncol(W), ncol = ncol(X))
            upper_beta    = matrix( Inf, nrow = ncol(W), ncol = ncol(X))
            if (!is.null(C)) {
              start_gamma_0 =
                start_gamma = rep(   0, ncol(C))
              lower_gamma   = rep(-Inf, ncol(C))
              upper_gamma   = rep( Inf, ncol(C))
            }

            if (is.null(C)) {
              svdd =
                svd(attr(
                  mu(X, W, oneXotimesW, start_alpha, start_beta, 0),
                  "gradient")[, seq(1, ncol(W) * ncol(X)) + ncol(W)])$d

            } else {
              svdd =
                svd(attr(
                  mu(X, W, C, oneXotimesW, start_alpha, start_beta, start_gamma, 0),
                  "gradient")[, seq(1, ncol(W) * ncol(X)) + ncol(W)])$d
            }
            # Decreasing order improves convergence.
            # start_* is taken from preceding nls run
            # If failed, tried again with start_*_0
            sqrtlambdalist = c(
              exp(seq(log(max(svdd)) + 1,
                      log(min(svdd)) - 1,
                      length.out = 20)),
              0)

            for (sqrtlambda in sqrtlambdalist) {
              print(sqrtlambda)
              if (is.null(C)) {
                mod = nls(y ~ mu(X, W, oneXotimesW, alpha, beta, sqrtlambda),
                          data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                      X = as.matrix(X),
                                      W = W,
                                      oneXotimesW = oneXotimesW),
                          start = list(alpha = start_alpha,
                                       beta  = start_beta),
                          lower = c(lower_alpha,
                                    lower_beta),
                          upper = c(upper_alpha,
                                    upper_beta),
                          algorithm = "port",
                          control = nls.control(warnOnly = TRUE))
                if (! mod$convInfo$isConv) {
                  mod = nls(y ~ mu(X, W, oneXotimesW, alpha, beta, sqrtlambda),
                            data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                        X = as.matrix(X),
                                        W = W,
                                        oneXotimesW = oneXotimesW),
                            start = list(alpha = start_alpha_0,
                                         beta  = start_beta_0),
                            lower = c(lower_alpha,
                                      lower_beta),
                            upper = c(upper_alpha,
                                      upper_beta),
                            algorithm = "port",
                            control = nls.control(warnOnly = TRUE))
                }
                if (! mod$convInfo$isConv) {
                  next()
                }
                start_alpha = coef(mod)[seq(1, ncol(W))]
                start_beta  = matrix(coef(mod)[seq(1, ncol(W) * ncol(X)) + ncol(W)],
                                     nrow = ncol(W), ncol = ncol(X))
              } else {
                mod = nls(y ~ mu(X, W, C, oneXotimesW, alpha, beta, gamma, sqrtlambda),
                          data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                      X = as.matrix(X),
                                      W = W,
                                      C = as.matrix(C),
                                      oneXotimesW = oneXotimesW),
                          start = list(alpha = start_alpha,
                                       beta  = start_beta,
                                       gamma = start_gamma),
                          lower = c(lower_alpha,
                                    lower_beta,
                                    lower_gamma),
                          upper = c(upper_alpha,
                                    upper_beta,
                                    upper_gamma),
                          algorithm = "port",
                          trace = TRUE,
                          control = nls.control(warnOnly = TRUE))
                if (mod$convInfo$isConv) {
                  start_alpha = coef(mod)[seq(1, ncol(W))]
                  start_beta  = matrix(coef(mod)[seq(1, ncol(W) * ncol(X)) + ncol(W)],
                                       nrow = ncol(W), ncol = ncol(X))
                  start_gamma = coef(mod)[seq(1, ncol(C)) + ncol(W) * (ncol(X) + 1)]
                } else {
                  start_alpha = coef(mod)[seq(1, ncol(W))]
                  mod = nls(y ~ mu(X, W, C, oneXotimesW, start_alpha, beta, gamma, sqrtlambda,
                                   gradientwithoutalpha = TRUE),
                            data = list(y = c(y, rep(0, ncol(X) * ncol(W))),
                                        X = as.matrix(X),
                                        W = W,
                                        C = as.matrix(C),
                                        oneXotimesW = oneXotimesW),
                            start = list(beta  = start_beta,
                                         gamma = start_gamma),
                            algorithm = "port",
                            trace = TRUE,
                            control = nls.control(warnOnly = TRUE))
                  if (! mod$convInfo$isConv) {
                    next()
                  } else {
                    start_beta  = matrix(coef(mod)[seq(1, ncol(W) * ncol(X))],
                                         nrow = ncol(W), ncol = ncol(X))
                    start_gamma = coef(mod)[seq(1, ncol(C)) + ncol(W) * ncol(X)]
                  }
                }
              }
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
              P =
                x[1:length(y), ] %*%
                solve(t(x) %*% x) %*%
                t(x[1:length(y), ])
              dof = sum(diag(P))
              RSS = sum((residuals(mod)[1:length(y)])^2)
              if (sqrtlambda == 0) {
                sigma2 = RSS / (length(y) - dof)
              }
              GCV = length(y) * RSS / (length(y) - dof)^2
              # AIC = 2 * dof +
              #   RSS / sigma2 +
              #   length(y) * log(2 * pi * sigma2)
              # BIC = log(length(y)) * dof +
              #   RSS / sigma2 +
              #   length(y) * log(2 * pi * sigma2)
              GCVstats = rbind(
                GCVstats,
                data.frame(sqrtlambda = sqrtlambda,
                           dof = dof,
                           GCV = GCV))
              start_alphas = c(start_alphas, list(start_alpha))
              start_betas  = c(start_betas,  list(start_beta))
              start_gammas = c(start_gammas, list(start_gamma))
            }
            if (0 %in% GCVstats$sqrtlambda) {
              i = which.min(GCVstats$GCV)
              sqrtlambda = GCVstats$sqrtlambda[i]
              dof = GCVstats$dof[i]
              start_alpha = start_alphas[[i]]
              start_beta  = start_betas[[i]]
              start_gamma = start_gammas[[i]]
              res = data.frame(estimate = c(start_alpha, start_beta, start_gamma))

              # Wald test
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
              sigma2Hstar =
                t(x[1:length(y), ]) %*%
                x[1:length(y), ]
              sigma2Hstarlambdainv = solve(t(x) %*% x)
              SE = sqrt(sigma2 * diag(sigma2Hstarlambdainv %*%
                                        sigma2Hstar %*%
                                        sigma2Hstarlambdainv))
              res$statistic = res$estimate / SE
              res$p.value = pt(- abs(res$statistic), df = length(y) - dof) * 2
              res$celltypeterm = c(colnames(oneXotimesW), colnames(C))
            } else {
              res =
                data.frame(estimate     = NA,
                           statistic    = NA,
                           p.value      = NA,
                           celltypeterm = c(colnames(oneXotimesW), colnames(C)))

            }
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
