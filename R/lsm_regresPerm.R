#' @title Regression tests for symptom mapping (permutation p-vlaues)
#'
#' @description
#' Lesion to symptom mapping performed on a prepared matrix.
#' Regressions are performed between behavior and each column in the
#' lesmat matrix. This function relies on the lmPerm package
#' to run. The number of permutations required to reach
#' stable p-values is estalbished automatically. For this reason,
#' the user cannot specify a predefined number of permutations.
#'
#' @param lesmat matrix of voxels (columns) and subjects (rows).
#' @param behavior vector of behavioral scores.
#'
#' @return List with vectors of statistic, pvalue, and zscore.
#'
#' @author Dorian Pustina
#'
#' @export
#' @importFrom lmPerm lmp

lsm_regresPerm <- function(lesmat, behavior) {
  statistic = pvalue = rep(NA, ncol(lesmat))

  if (! 'lmPerm' %in% rownames(installed.packages())) stop('regresPerm requires lmPerm package. Try installing with install.packages("lmPerm")')

  temp = lmPerm::lmp(lesmat ~ behavior, perm = 'Prob')
  temp = ANTsRCore::bigLMStats(temp)
  statistic = temp$beta.t
  pvalue = temp$beta.pval
  zscore = qnorm(pvalue, lower.tail = FALSE)

  return(list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore))
}
