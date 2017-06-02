#' lsm_BMfast
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
#' @return List with vectors of statistic, pvalue, zscore.
#' If FWER is selected, perm.vector and perm.FWERthresh
#' will also be returned.
#'
#' @author Dorian Pustina
#'
#' @export
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



  p.value = rep(1, length(statistic))
  zscore =  rep(0, length(statistic))

  if (!statOnly) {
    if ((alternative == "less") | (alternative == "l")) {
      p.value = pt(statistic, dof)
    }
    else if ((alternative == "greater") | (alternative == "g")) {
      p.value = pt(statistic, dof, lower.tail=F)
    }
    else {
      alternative = "two.sided"
      p.value = 2 * pt(abs(statistic), dof, lower.tail=F)
        # apply( rbind(   pt(statistic, dof, lower.tail=F)  , pt(statistic, dof, lower.tail=T) ), 2, min)
    }

    # compute zscores
    zscore = qnorm(dt(statistic,dof))
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
      if (alternative == 'greater') maxvec[p] = tempperm[v]
      if (alternative == 'less') maxvec[p] = tempperm[ length(tempperm) - v + 1 ]
      if (alternative == 'two.sided') maxvec[p] = abs(tempperm)[v]
    }

    # threshold statistic, abs() needed to consider when peak statistic is negative
    if (alternative == 'greater') FWEthresh = quantile(maxvec, probs = (1-pThreshold))
    if (alternative == 'less') FWEthresh = quantile(maxvec, probs = (pThreshold))
    if (alternative == 'two.sided') FWEthresh = quantile(maxvec, probs = (1-pThreshold)) # no need to take half pThreshold because we recorded the absolute, all top values are in the high side anyway

    if (alternative == 'greater') statistic[statistic < FWEthresh] = 0
    if (alternative == 'less') statistic[statistic > FWEthresh] = 0
    if (alternative == 'less') statistic[abs(statistic) > FWEthresh] = 0
  }



  # return outcome
  output = list(statistic=statistic)
  if (!statOnly) {
    output$pvalue = p.value
    output$zscore = zscore
  }

  if (FWERperm) output$perm.FWEthresh = FWEthresh
  if (FWERperm) output$perm.vector = maxvec

  return(output)
}
