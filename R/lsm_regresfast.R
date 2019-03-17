#' lsm_regresfast
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' Regressions are performed between behavior and each column
#' of the lesmat matrix. Fast function based on compiled code.
#'
#' @param lesmat matrix of voxels (columns) and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param covariates (default=NA) vector of matrix of covariates.
#' @param FWERperm logical (default=FALSE) whether to run permutation
#' based FWER thresholding.
#' @param v (default=1) what voxel to record for FWER thresholding.
#' @param clusterPerm logical (default=FALSE), whether to perform
#' permutation based cluster thresholding.
#' @param mask (default=NA) antsImage reference mask used for
#' cluster computations.
#' @param voxindx (default=NA) indices of voxels to put in mask
#' @param samplemask (default=NA) antsImage used to extract voxels
#' back in a matrix.
#' @param pThreshold (default=0.05) Voxel-wise threshold.
#' @param clusterPermThreshold (default=0.05) threshold for cluster
#' selection after obtaining cluster size distrubution.
#' @param nperm Number of permutations to perform when needed.
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
#' \item\code{perm.clusterThreshold} - (optional) permutation threshold established
#' from the distribution of \code{perm.vector}
#'
#' }
#'
#' @examples{
#' set.seed(123)
#' lesmat = matrix(rbinom(200,1,0.5), ncol=2)
#' set.seed(123)
#' behavior = rnorm(100)
#' result = lsm_regresfast(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
#' @export


lsm_regresfast <- function(lesmat, behavior, covariates=NA,
                           FWERperm=F, nperm=1000, v=1, pThreshold=0.05,
                           clusterPerm = F, mask = NA, voxindx = NA, samplemask= NA, clusterPermThreshold=0.05,
                           showInfo=T, ...) {

#   if (!usePkg('Rcpp') | !usePkg('RcppArmadillo'))
#     stop('Packages Rcpp and RcppArmadillo must be installed to use BMfast.')
#   # below source still doesn't work in R CMD call
#   if (length(find('regresfast', mode='function')) == 0) {
#     Rcpp::sourceCpp(file.path( getSrcDirectory(lsm_BMfast) , 'regresfast.cpp'))
#   }

  # check covariates are ok
  hascovariate = F
  if (!is.na(covariates[1])) {
    if (class(covariates) != 'matrix') covariates=as.matrix(covariates)
    if (nrow(covariates) != nrow(lesmat)) stop('Mismatch rows in covariates and lesion matrix.')
    hascovariate = T
    ncovariates = ncol(covariates)
  } else {
    covariates = as.matrix(rep(0, nrow(lesmat)))
    ncovariates = 0
  }


  # run default analysis
  temp = regresfast(lesmat, behavior, covariates, hascovariate)
  statistic = temp$statistic
  pvalue = pt( abs(statistic) * -1, temp$n - temp$kxmat, lower.tail = T, log.p = F) * 2


  # user has requested permutation thresholding
  if (FWERperm | clusterPerm) {

    # prepare coviariate for Freedman-Lane approach of permutations
    if (hascovariate) {
      temp = lm(behavior ~ covariates)
      covar.coef = temp$coefficients
      covar.resid = residuals(temp)
    }

    # run permutations
    maxvec = rep(NA, nperm)
    for (p in 1:nperm) {

      # keep time for estimating the runtime
      if (p == 1) tic = Sys.time()
      if (p == 2 & showInfo) {
        toc = tic + (nperm * onerun)
        expect = paste(round(as.double(difftime(toc,tic)),1), units(difftime(toc,tic)))
        covarinfo = ''
        if (hascovariate) covarinfo = 'with covariates'
        mytype = ifelse(FWERperm, 'FWERperm', 'clusterPerm')
        printInfo(paste('\n       ',mytype, nperm, 'permutations', covarinfo,' - expected run =', expect), type='middle')
      }

      # permute behavior
      rindx = sample(1:length(behavior), replace = F)
      if (hascovariate) { # FREEDMAN-LANE
        behavior.permuted = covariates%*%covar.coef[-1] + covar.coef[1] + covar.resid[rindx]
      } else { # NO COVARIATES
        behavior.permuted = behavior[rindx]
      }

      # run this permutation
      tempperm = regresfast(lesmat, behavior.permuted, covariates, hascovariate)

      # record this value
      if (FWERperm) {

        maxvec[p] = sort(abs(tempperm$statistic), decreasing = T, method='quick')[v]

      } else if (clusterPerm) {

        if (any(is.na( c(mask, voxindx, samplemask, clusterPermThreshold) ))) stop('Missing proper inputs to run clusterPerm')

        # mask, voxindx, pThreshold, clusterPermThreshold, samplemask, needed from lesymap
        perm.pvalue = pt( abs(tempperm$statistic) * -1, tempperm$n - tempperm$kxmat) * 2
        tempperm$statistic[perm.pvalue > pThreshold] = 0
        permstat = makeImage(mask, abs(tempperm$statistic[voxindx]))
        lstat = labelStats( mask, labelClusters(permstat, minClusterSize = 1,
                                                   minThresh = 0.01, maxThresh = Inf, fullyConnected = T))
        if (nrow(lstat)>1) maxvec[p] = max(lstat$Volume[lstat$LabelValue!=0])
        else maxvec[p] = 0
        invisible(gc()) # this steals time but might be necessary in linux machines due to memory leaks with antsImages
      }

      # keep time for estimating runtime
      if (p == 1) onerun = as.double(difftime(Sys.time(),tic, units = 'sec'))
    }

    # permutation done, threshold statistics
    if (FWERperm) {

      # threshold statistic from permutation distribution, abs() needed to consider when peak statistic is negative
      FWEthresh = quantile(maxvec, probs = (1-pThreshold))
      statistic[abs(statistic) < FWEthresh] = 0
      pvalue = pt( abs(statistic) * -1, temp$n - temp$kxmat, lower.tail = T, log.p = F) * 2

    } else if (clusterPerm) {

      pvalue = pt( abs(statistic) * -1, temp$n - temp$kxmat, lower.tail = T, log.p = F) * 2
      statistic[pvalue > pThreshold] = 0 # threshold statistic
      # real statistic image
      realstat = makeImage(mask, abs(statistic[voxindx]) ) # create temporary stat image
      clust.thresh = quantile(maxvec, probs = (1-clusterPermThreshold))
      realclust = labelClusters(realstat, minClusterSize = as.integer( round(clust.thresh,0) ),
                                minThresh = 0.01, maxThresh = Inf, fullyConnected = T)
      clustmask = thresholdImage(realclust, 0.1, Inf) # binarize surviving clusters
      clustbin = imageListToMatrix(list(clustmask), samplemask)[1,] # get array of surviving voxels
      statistic = statistic*clustbin # eliminate voxels from other clusters
      #       realstat = realstat*clustmask
      #       statistic = imageListToMatrix(list(realstat), samplemask)[1,]
      pvalue = pt( abs(statistic) * -1, temp$n - temp$kxmat, lower.tail = T, log.p = F) * 2

    }


  }

  zscore = qnorm(pvalue, lower.tail = FALSE)
  # zscore = qt(pvalue, length(behavior) - 2 - ncovariates)

  # return outcome
  output = list(statistic=statistic,
                pvalue=pvalue,
                zscore=zscore)
  if (FWERperm | clusterPerm) output$perm.vector = maxvec
  if (FWERperm) output$perm.FWERthresh = FWEthresh
  if (clusterPerm) output$perm.clusterThreshold = as.integer( round(clust.thresh,0) )

  return(output)
}
