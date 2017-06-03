#' lsm_BM
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' Brunner-Munzel tests are performed using each column
#' of the matrix to split the behavioral scores in two
#' groups.
#'
#' @param lesmat binary matrix (0/1) of voxels (columns)
#' and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param permuteNthreshold (default=9) Voxels lesioned in less than
#' this number will undergo permutation based thresholding.
#' See Medina et al 2010.
#' @param nperm Number of permutations to perform when needed.
#' @param alternative (default="greater") It is assumed that
#' healthy voxels (0) have greater behavioral scores. If your
#' data follow an inverted relationship choose "less" or
#' "two.sided".
#' @param showInfo display info messagges when running the function.
#' @param ... other arguments received from \code{\link{lesymap}}.
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
lsm_BM <- function(lesmat, behavior, permuteNthreshold=9, nperm=10000, alternative='greater', showInfo=T, ...) {

  statistic = pvalue = rep(NA, ncol(lesmat))

  # find voxels that need permutation
  lesmat.colsums = colSums(lesmat)
  permindx = (lesmat.colsums)<=permuteNthreshold
  permindx = (nrow(lesmat) - lesmat.colsums)<=permuteNthreshold | permindx

  # run standard non permuted Brunner Munzel
  # need to run on all voxels to get degrees of freedom for zscoring
  temp = BM(lesmat, behavior, alternative=alternative)
  statistic[!permindx] = temp$statistic[!permindx]
  dof = temp$dof

  # fix problem of infinite dof
  dof[dof==Inf] = .Machine$double.xmax


  if ((alternative == "less") | (alternative == "l")) {
    pvalue = pt(statistic, dof)
  }
  else if ((alternative == "greater") | (alternative == "g")) {
    pvalue = pt(statistic, dof, lower.tail=F)
  }
  else {
    alternative = "two.sided"
    pvalue = 2 * pt(abs(statistic), dof, lower.tail=F)
  }


  # run permutation test on selected voxels
  if (sum(permindx) > 0) { # run only if any voxel needs permutated brunner-munzel

    if (showInfo) cat(paste0('\n        running permutation on ', sum(permindx),' voxels below permuteNthreshold' ))
    if (! 'nparcomp' %in% rownames(installed.packages())) stop('Permutation not possible without the nparcomp package. Try installing with install.packages("nparcomp")')

    output = apply(lesmat[,permindx], 2, function(x) {
      temp = nparcomp::npar.t.test(behavior~group, alternative=alternative,
                                   data=data.frame(behavior=behavior, group=x),
                                   method='permu', nperm=nperm, info=F)
      # careful, nparcomp has correct alternative pvalues, but inverted statistics
      return(list(
        stat= -temp$Analysis$Statistic[1],
        pval= temp$Analysis$p.value[1]))
    }
    )

    temp = unlist(output)
    statistic[permindx] = temp[seq(1,length(temp),by=2)]
    pvalue[permindx] = temp[seq(2,length(temp),by=2)]
    rm(temp, output)
  } else {
    if (showInfo) cat(paste0('\n        No permutation needed, all voxels above permuteNthreshold.' ))
  }









  # compute zscores
  zscore = qnorm(pvalue)
  # use below if you get values -Inf in zscores from p-values
  # which happens because of precision limitations
  # zscore = qnorm(dt(statistic,dof))

  return(list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore))
}
