#' Small Subset of GSE42861 Dataset From GEO
#'
#' The dataset includes 336 rheumatoid arthritis cases and 322 controls.
#' A subset of 500 CpG sites were randomly selected from the original EWAS dataset.
#'
#' @docType data
#'
#' @usage data(GSE42861small)
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861}{GEO}
#'
#' @seealso ctassoc
#' @examples
#' data(GSE42861small)
#' X = GSE42861small$X
#' W = GSE42861small$W
#' Y = GSE42861small$Y
#' Y = Y[seq(1, 20), ] # for brevity
#' C = GSE42861small$C
#' result = ctassoc(X, W, Y, C = C)
#' result$coefficients
"GSE42861small"
