#' lsm_regres
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' Regressions are performed between behavior and each column in the
#' lesmat matrix.
#'
#' @param lesmat matrix of voxels (columns) and subjects (raws).
#' @param behavior vector of behavioral scores.
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{statistic} - vector of statistical values
#' \item\code{pvalue} - vector of pvalues
#' \item\code{zscore} - vector of zscores
#' }
#'
#' @author Dorian Pustina
#'
#' @export
lsm_regres <- function(lesmat, behavior) {
  statistic = pvalue = rep(NA, ncol(lesmat))

  temp = lm(lesmat ~ behavior)
  temp = bigLMStats(temp)
  statistic = temp$beta.t
  pvalue = temp$beta.pval
  zscore = qnorm(pvalue)

  return(list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore))
}
