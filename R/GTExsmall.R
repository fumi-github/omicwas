#' Small Subset of GTEx Dataset
#'
#' The dataset includes gene expression measured in whole blood for 389 samples.
#' A subset of 500 genes were randomly selected from the original dataset.
#'
#' @docType data
#'
#' @usage data(GTExsmall)
#'
#' @source \href{https://gtexportal.org}{GTEx}
#'
#' @seealso ctassoc
#' @examples
#' data(GTExsmall)
#' X = GTExsmall$X
#' Y = GTExsmall$Y + 1
#' Y = Y[seq(1, 20), ] # for brevity
#' W = GTExsmall$W
#' C = GTExsmall$C
#' result = ctassoc(X, W, Y, C = C)
#' result$coefficients
"GTExsmall"
