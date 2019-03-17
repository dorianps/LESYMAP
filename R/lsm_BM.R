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
#'  \item\code{statistic} - vector of statistical values
#'  \item\code{pvalue} - vector of pvalues
#'  \item\code{zscore} - vector of zscores
#' }
#'
#' @examples{
#' set.seed(123)
#' lesmat = matrix(rbinom(200,1,0.5), ncol=2)
#' set.seed(123)
#' behavior = rnorm(100)
#' result = lsm_BM(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
#' @export
# @importFrom nparcomp npar.t.test
lsm_BM <- function(lesmat, behavior, permuteNthreshold=9, nperm=10000,
                   alternative='greater', showInfo=TRUE, ...) {

  if (showInfo) warning('The BM method is deprecated and will be removed from future LESMAP versions. Please use method=\'BMfast\'.')
  statistic = pvalue = rep(NA, ncol(lesmat))

  # find voxels that need permutation
  lesmat.colsums = colSums(lesmat)
  permindx = (lesmat.colsums)<=permuteNthreshold
  permindx = (nrow(lesmat) - lesmat.colsums)<=permuteNthreshold | permindx

  # run standard non permuted Brunner Munzel
  # need to run on all voxels to get degrees of freedom for zscoring
  # temp = BM(lesmat, behavior, alternative=alternative)

  temp = BM(lesmat, behavior)
  statistic[!permindx] = temp$statistic[!permindx]
  dof = temp$dof

  # fix problem of infinite dof
  dof[dof==Inf] = .Machine$double.xmax


  if ((alternative == "less") | (alternative == "l")) {
    pvalue = pt(statistic, dof, lower.tail=TRUE)
  } else if ((alternative == "greater") | (alternative == "g")) {
    pvalue = pt(statistic, dof, lower.tail=FALSE)
  } else {
    alternative = "two.sided"
    pvalue = 2 * pt(abs(statistic), dof, lower.tail=FALSE)
  }


  # run permutation test on selected voxels
  if (sum(permindx) > 0) { # run only if any voxel needs permutated brunner-munzel

    if (showInfo) printInfo(paste0('\n        running ', nperm, ' permutations on ', sum(permindx),' voxels below permuteNthreshold' ), type='middle')
    if (nperm < 20000) warning('Number of permutations too small, consider increasing it.')
    if (! 'nparcomp' %in% rownames(installed.packages())) {
      stop('Permutation not possible without the nparcomp package. Try installing with install.packages("nparcomp")')
    }

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

    # fixing pvalue == 0, must be a problem in nparcomp
    pvalue[pvalue == 0] = 1/(nperm+1)

    rm(temp, output)
  } else {
    if (showInfo) printInfo(paste0('\n        No permutation needed, all voxels above permuteNthreshold.' ), type='middle')
  }


  # compute zscores after permutation has updated pvalues
  if ((alternative == "less") | (alternative == "l")) {
    zscore = qnorm(pvalue, lower.tail=TRUE)
  }
  else if ((alternative == "greater") | (alternative == "g")) {
    zscore = qnorm(pvalue, lower.tail=FALSE)
  } else {
    neg = statistic<0
    pos = statistic>0
    zscore = statistic*0
    zscore[neg] = qnorm(pvalue[neg], lower.tail=TRUE)
    zscore[pos] = qnorm(pvalue[pos], lower.tail=FALSE)
  }

  # useless commands/comments
  # compute zscores
  # zscore = qnorm(pvalue)
  # zscore = qt(pvalue, dof)
  # use below if you get values -Inf in zscores from p-values
  # which happens because of precision limitations
  # zscore = qnorm(dt(statistic,dof))

  # trying to avoid infitive values
  zscore[is.nan(zscore)] = .Machine$double.xmax
  zscore[is.na(zscore)] = .Machine$double.xmax
  zscore[zscore==Inf] = .Machine$double.xmax
  zscore[zscore==-Inf] = -.Machine$double.xmax


  return(list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore))
}
