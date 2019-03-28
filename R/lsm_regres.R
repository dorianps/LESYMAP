#' @title Regression tests for symptom mapping (slow)
#'
#' @description
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
#' @examples{
#' set.seed(123)
#' lesmat = matrix(rbinom(200,1,0.5), ncol=2)
#' set.seed(123)
#' behavior = rnorm(100)
#' result = lsm_regres(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
#' @export
lsm_regres <- function(lesmat, behavior) {
  statistic = pvalue = rep(NA, ncol(lesmat))

  temp = lm(lesmat ~ behavior)
  dof = temp$df.residual
  temp = bigLMStats(temp)
  statistic = temp$beta.t
  pvalue = temp$beta.pval
  #zscore = qt(pvalue, length(behavior)-2)
  zscore = qnorm(pvalue, lower.tail = FALSE)

  return(list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore))
}
