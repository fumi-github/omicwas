##' Fitting reduced-rank ridge regression with given rank and shrinkage penalty
##'
##' Fitting reduced-rank ridge regression with given rank and shrinkage penalty
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param nrank an integer specifying the desired rank
##' @param lambda tunging parameter for the ridge penalty
##' @param coefSVD logical indicating the need for SVD for the
##'   coeffient matrix int the output
##' @return S3 \code{rrr} object, a list consisting of
##'   \item{coef}{coefficient of rrs}
##'   \item{coef.ls}{coefficient of least square}
##'   \item{fitted}{fitted value of rrs}
##'   \item{fitted.ls}{fitted value of least square}
##'   \item{A}{right singular matrix}
##'   \item{Ad}{sigular value vector}
##'   \item{nrank}{rank of the fitted rrr}
##' @examples
##' library(rrpack)
##' Y <- matrix(rnorm(400), 100, 4)
##' X <- matrix(rnorm(800), 100, 8)
##' rfit <- rrs.fit(Y, X)
##' @references
##'
##' Mukherjee, A. and Zhu, J. (2011) Reduced rank ridge regression and its
##' kernal extensions.
##'
##' Mukherjee, A., Chen, K., Wang, N. and Zhu, J. (2015) On the degrees of
##' freedom of reduced-rank estimators in multivariate
##' regression. \emph{Biometrika}, 102, 457--477.
##'
##' @importFrom MASS ginv
##' @export
rrs.fit <- function(Y,
                    X,
                    nrank = min(ncol(Y), ncol(X)),
                    lambda = 1,
                    coefSVD = FALSE)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  S_yx <- t(Y) %*% X
  ## This is a key difference
  S_xx <- t(X) %*% X + lambda * diag(p)

  ## S_xx_inv <- tryCatch(solve(S_xx+0.1*diag(p)),error=function(e)ginv(S_xx))
  ## S_xx_inv <- ginv(S_xx)
  ## Use the Woodbury matrix identity
  if (lambda != 0) {
    S_xx_inv <- 1 / lambda * diag(p) -
      lambda ^ (-2) * t(X) %*% ginv(diag(n) + lambda ^ (-1) * X %*%
                                      t(X)) %*% X
  } else{
    S_xx_inv <- ginv(S_xx)
    if (sum(is.na(S_xx_inv)) > 0) {
      S_xx_inv <- solve(S_xx + 0.1 * diag(p))
    }
  }

  C_ls <- S_xx_inv %*% t(S_yx)

  ypy.svd <- TRUE
  ##if(ypy.svd){
  ##This is another key difference
  XC <- rbind(X, sqrt(lambda) * diag(p)) %*% C_ls
  svdXC <- tryCatch(
    svd(XC, nu = nrank, nv = nrank),
    error = function(e)
      2)
  if (tryCatch(
    svdXC == 2,
    error = function(e)
      3) == 3) {
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
  } else{
    ypy.svd <- FALSE
  }
  #}
  if (!ypy.svd) {
    SS <- S_yx %*% C_ls
    SS <- (SS + t(SS)) / 2
    eigenSS <- eigen(SS, symmetric = TRUE)
    A <- as.matrix(eigenSS$vectors[, 1:nrank])
    Ad <- eigenSS$values[1:nrank]
  }

  AA <- A %*% t(A)
  C_rr <- C_ls %*% AA

  ##    if(c.svd){
  ##      svd_C <- svd(C_rr,nv=nrank,nu=nrank)
  ##      U <- as.matrix(svd_C$u[,1:nrank])
  ##      V <- as.matrix(svd_C$v[,1:nrank])
  ##      D <- diag(svd_C$d[1:nrank],nrow=nrank)
  ##
  ##      ####return ls estimator C_ls, reduced-rank estimator C_rr
  ##      ####return SVD of C_rr
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,U=U,V=V,D=D,C=C_rr,rank=nrank)
  ##    }else{
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,C=C_rr,rank=nrank)
  ##    }

  ret <- list(
    call = Call,
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    nrank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  class(ret) <- "rrs.fit"
  ret
}
