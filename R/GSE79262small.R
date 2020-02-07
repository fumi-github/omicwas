#' Small Subset of GSE79262 Dataset From GEO
#'
#' The dataset includes 53 samples.
#' A subset of 737 CpG sites and 3624 SNPs within Chr1:100,000,000-110,000,000
#' were selected from the original EWAS dataset.
#' DNA methylation was measured in T cells.
#' The estimated proportion of CD4T, CD8T, NK cells are saved in W.
#'
#' @docType data
#'
#' @usage data(GSE79262small)
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79262}{GEO}
#'
#' @seealso ctcisQTL
#' @examples
#' data(GSE79262small)
#' X    = GSE79262small$X
#' Xpos = GSE79262small$Xpos
#' W    = GSE79262small$W
#' Y    = GSE79262small$Y
#' Ypos = GSE79262small$Ypos
#' C    = GSE79262small$C
#' X    = X[seq(1, 3001, 100), ] # for brevity
#' Xpos = Xpos[seq(1, 3001, 100)]
#' Y    = Y[seq(1, 501, 100), ]
#' Ypos = Ypos[seq(1, 501, 100)]
#' ctcisQTL(X, Xpos, W, Y, Ypos, C = C)
"GSE79262small"
