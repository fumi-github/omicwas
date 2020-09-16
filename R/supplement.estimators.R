# Supplement new estimators that can be computed by available estimators
# and return summary matrix.
# new_estimators = estimatormatrix %*% available_estimators
.supplement.estimators = function (mlm, estimatormatrix) {
  df.residual = mlm$df.residual
  qr = mlm$qr
  rank = mlm$rank
  weights = NULL
  result =
    purrr::pmap(
      list(
        coefficients = plyr::alply(mlm$coefficients, 2),
        residuals = plyr::alply(mlm$residuals, 2)),
      function(coefficients, residuals) {
        x = list(
          df.residual = df.residual,
          qr = qr,
          rank = rank,
          weights = weights,
          coefficients = coefficients,
          residuals = residuals)
        v = .Vcov.lm(x)
        estimate = c(coefficients,
                     estimatormatrix %*% coefficients)
        SE = sqrt(c(diag(v),
                    diag(estimatormatrix %*% v %*% t(estimatormatrix))))
        statistic = estimate / SE
        p.value = pt(abs(statistic), df = df.residual, lower.tail = FALSE) *2
        return(data.frame(
          celltypeterm = c(colnames(estimatormatrix), rownames(estimatormatrix)),
          estimate = estimate,
          statistic = statistic,
          p.value = p.value))
      })
  names(result) = colnames(mlm$coefficients)
  return(result)
}

# Standard vcov function returns a large matrix for mlm object.
# Instead, apply the following to each lm object.
# Copied from vcov.R
.Vcov.lm = function(object, ...) {
  if (p <- object$rank) {
    p1 = seq_len(p)
    rss = if (is.null(w <- object$weights)) {
      sum(object$residuals^2)
    } else {
      sum(w * object$residuals^2)
    }
    covmat = rss * chol2inv(object$qr$qr[p1, p1, drop = FALSE])/
      object$df.residual
    nm = names(object$coefficients)
    dimnames(covmat) = list(nm, nm)
    return(covmat)
  } else return(numeric(0))
}
