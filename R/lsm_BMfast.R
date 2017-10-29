#' lsm_BMfast
#'
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' Brunner-Munzel tests are performed using each column
#' of the matrix to split the behavioral scores in two
#' groups. This function relies on a compiled version for
#' fast processing.
#'
#' @param lesmat binary matrix (0/1) of voxels (columns)
#' and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param permuteNthreshold (default=9) Voxels lesioned in less than
#' this number will undergo permutation based thresholding.
#' See Medina et al 2010.
#' @param nperm (default=1000) Number of permutations to perform when needed.
#' @param alternative (default="greater") It is assumed that
#' healthy voxels (0) have greater behavioral scores. If your
#' data follow an inverted relationship choose "less" or
#' "two.sided".
#' @param statOnly logical (default=FALSE), skips some computations,
#' don't use unless you know it's effects
#' @param FWERperm logical (default=FALSE) whether to perform permutation
#' based FWER thresholding.
#' @param v (default=1) which voxel to record at each permutation, first
#' or other voxels (i.e., v=10 for 10 highest voxel)
#' @param pThreshold (default=0.05) what threshold to use for FWER
#' @param showInfo display info messagges when running the function.
#' @param ... other arguments received from \code{\link{lesymap}}.
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{statistic} - vector of statistical values
#' \item\code{pvalue} - vector of pvalues
#' \item\code{zscore} - vector of zscores
#' \item\code{perm.vector} - (optional) vector of permuted statistics
#' \item\code{perm.FWERthresh} - (optional) permutation threshold established
#' from the distribution of \code{perm.vector}
#' }
#' @export
#' @examples{
#' set.seed(123)
#' lesmat = matrix(rbinom(200,1,0.5), ncol=2)
#' set.seed(123)
#' behavior = rnorm(100)
#' result = lsm_BMfast(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
lsm_BMfast <- function(lesmat, behavior, permuteNthreshold=9, alternative="greater",
                       statOnly = F, nperm=1000,
                       FWERperm=F, v=1, pThreshold=0.05,
                       showInfo=F,...) {


  # check the assumption of min subjects in BM is not violated
  # find voxels that need permutation
  lesmat.colsums = colSums(lesmat)
  permindx = (lesmat.colsums)<=permuteNthreshold
  permindx = (nrow(lesmat) - lesmat.colsums)<=permuteNthreshold | permindx
  if (sum(permindx) > 0 & !FWERperm) stop(paste(sum(permindx), 'voxels need BM permutation. This is not possible in BMfast, please choose regular BM.'))

  if (FWERperm) statOnly = T

  # main BM run
  tic = Sys.time()
  temp = BMfast2(lesmat, behavior, computeDOF = !FWERperm)
  onerun = as.double(difftime(Sys.time(),tic, units = 'sec'))
  statistic = temp$statistic
  dof = temp$dfbm

  # fix problem of infinite dof
  dof[dof==Inf] = .Machine$double.xmax



  pvalue = rep(1, length(statistic))
  zscore =  rep(0, length(statistic))

  if (!statOnly) {
    if ((alternative == "less") | (alternative == "l")) {
      pvalue = pt(statistic, dof, lower.tail=TRUE)
      zscore = qnorm(pvalue, lower.tail=TRUE)
      #zscore = qt(pvalue, dof, , lower.tail=TRUE)
    }
    else if ((alternative == "greater") | (alternative == "g")) {
      pvalue = pt(statistic, dof, lower.tail=FALSE)
      zscore = qnorm(pvalue, lower.tail=FALSE)
      #zscore = qt(pvalue, dof, lower.tail=FALSE)
    }
    else {
      alternative = "two.sided"
      pvalue = 2 * pt(abs(statistic), dof, lower.tail=FALSE)
      zscore = qnorm(pvalue, lower.tail=FALSE)
      #zscore = qt(pvalue, dof, , lower.tail=FALSE)
    }

    #' Note on zscores
    #' qnorm gives same values as MRIcron
    #' and relies on the normal distribution.
    #' however, we are computing t-scores, and
    #' should have relied on that distribution,
    #' which is the t-score itself.

  }

  if (FWERperm) {
    if (showInfo) {
      toc = tic + (nperm * onerun)
      expect = paste(round(as.double(difftime(toc,tic)),1), units(difftime(toc,tic)))
      cat(paste('\n       FWERperm', nperm, 'permutations, expected run =', expect))
    }
    # run permutations
    maxvec = rep(NA, nperm)
    for (p in 1:nperm) {
      tempperm = BMfast2(lesmat, sample(behavior), computeDOF=FALSE)
      tempperm = sort(tempperm$statistic, decreasing = TRUE, method='quick')
      # if we expect greater, get max value
      if (alternative == 'greater') maxvec[p] = tempperm[v]
      # if we expect less, get min value
      if (alternative == 'less') maxvec[p] = tempperm[ length(tempperm) - v + 1 ]
      # if we expect twosided, get most extreme value
      if (alternative == 'two.sided') {
        thismax = tempperm[v]
        thismin = tempperm[ length(tempperm) - v + 1 ]
        maxvec[p] = ifelse(abs(thismin) > abs(thismax), thismin, thismax)
      }
    }

    if (alternative == 'greater') FWEquantile = (1-pThreshold)
    if (alternative == 'less') FWEquantile = (pThreshold)
    if (alternative == 'two.sided') FWEquantile = c(pThreshold/2,(1-pThreshold/2))
        
    # threshold statistic, abs() needed to consider when peak statistic is negative
    if (alternative == 'greater') FWEthresh = quantile(maxvec, probs = FWEquantile)
    if (alternative == 'less') FWEthresh = quantile(maxvec, probs = FWEquantile)
    if (alternative == 'two.sided') FWEthresh = quantile(maxvec, probs = FWEquantile)

    if (alternative == 'greater') statistic[statistic < FWEthresh] = 0
    if (alternative == 'less') statistic[statistic > FWEthresh] = 0
    if (alternative == 'two.sided') statistic[statistic > FWEthresh[1] & statistic < FWEthresh[2] ] = 0
  }



  # return outcome
  output = list(statistic=statistic)
  output$pvalue = pvalue
  if (!statOnly) {
    output$zscore = zscore
  }

  if (FWERperm) output$perm.FWEquantile = paste(FWEquantile, collapse=' | ')
  if (FWERperm) output$perm.FWEthresh = paste(FWEthresh, collapse=' | ')
  if (FWERperm) output$perm.vector = maxvec
  

  return(output)
}
