#' Sparse canonical correlations for symptom mapping.
#'
#' Multivariate SCCAN adapted for lesion to symptom mapping purposes.
#' By default an optimization routine is used to find the best
#' \code{sparseness} value. If you specify sparseness manually, it
#' will be validated to find the cross-validated correlation that
#' can be obtained with that sparseness. You can skip the entire
#' optimization/validation by choosing \code{optimizeSparseness=FALSE}.
#' To understand SCCAN arguments, see \code{\link[ANTsR]{sparseDecom2}}.
#'
#' @param lesmat matrix of voxels (columns) and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param mask antsImage binary mask to put back voxels in image.
#' @param optimizeSparseness logical (default=TRUE) whether to
#' run the sparseness optimization routine. If FALSE, the default
#' sparseness value will be used. If sparseness is manually defined
#' this flag decides if cross validated correlations will be
#' computed for the defined sparseness.
#' @param validateSparseness logical (conditional default=TRUE) If
#' sparseness is manually defined, this flag decides if cross
#' validated correlations will be computed for the defined sparseness.
#' @param pThreshold (default=0.05) If cross validated
#' correlations show significance below this value
#' the results are considered null and an empty
#' map is returned.
#' @param showInfo logical (default=TRUE) display messages
#' @param sparseness (default=1) SCCAN parameter. Decides the proportion
#' of voxels that will receive a non-zero weight. A positive sparseness
#' will force the solution of each component to be one sided, i.e.,
#' voxels cannot have both positive and negative weights. A negative
#' sparseness allows dual sided solution, where some voxels can have
#' positive weights and other voxels can have negative weights. Setting
#' sparseness manually without running the optimization routing is not
#' recommended. For more, see \code{\link[ANTsR]{sparseDecom2}}.
#' @param sparseness.behav SCCAN parameter, what sparsness to use for
#' behavioral scores. Useful only if multiple behavioral scores are
#' passed. This argument is not optimized, you should not change it
#' if you are not familiar with SCCAN.
#' @param mycoption (default=1) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param robust (ddefault=1) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param nvecs (default=1) SCCAN parameter. Normally only
#' one eigenvector of weights is obtained in LESYMAP. Multiple
#' maps/eigenvectors can be retrieved for mapping full
#' deficit profiles in the future. For more, see
#' \code{\link[ANTsR]{sparseDecom2}}
#' @param cthresh (default=150) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param its (default=20) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param smooth (default=0.4) SCCAN parameter. Determines the
#' amount of smoothing of weights in image space performed by
#' \code{\link[ANTsR]{sparseDecom2}}. The current default value
#' is somewhat arbitrary, it was not determined through
#' systematic simulations.
#' @param npermsSCCAN (default=0) SCCAN permutations. In theory can be
#' used to determine if the cross-correlation between the two sides
#' (behavior and lesions) is not random. However, LESYMAP uses
#' k-fold validations, which are faster; this option has not been
#' tested. For more,  see \code{\link[ANTsR]{sparseDecom2}}.
#' @param maxBased (default=FALSE) SCCAN parameter. Removes voxels with
#' weights smaller than 10\% of the peak weight during internal SCCAN
#' iterations. Although similar to what is done in LESYMAP with standard
#' SCCAN results, this strategy follows a different route, and produces
#' different weights. The overall final result is, however, quite similar.
#' This method is faster then the standard SCCAN call in LESYMAP, but
#' has not been tested thoroughly. Note that the optimal sparseness
#' obtained with \code{maxBased=TRUE} is not optimal when switching to
#' \code{maxBased=FALSE}.
#' @param directionalSCCAN (default=TRUE) If TRUE, the upper and lower
#' bounds of sparseness search will be negative. A negative sparseness
#' permits positive and negative voxel weights, thus finding the
#' direction of the relationship with behavior.
#' @param ... other arguments received from \code{\link{lesymap}}.
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{statistic} - vector of statistical values
#' \item\code{pvalue} - vector of pvalues
#' \item\code{rawWeights.img} - image with raw SCCAN voxel weights
#' \item\code{sccan.eig2} - SCCAN weight(s) for behavior
#' column(s).
#' \item\code{sccan.ccasummary} - SCCAN summary of
#' projection correlations and permutation-derived pvalues
#' \item\code{optimalSparseness} - (if optimizeSparseness=TRUE) optimal
#' value found for sparseness
#' \item\code{CVcorrelation.stat} - (if optimizeSparseness=TRUE)
#' Correlation between true and predicted score with k-fold validation
#' using the optimal sparseness value
#' \item\code{CVcorrelation.pval} - (if optimizeSparseness=TRUE) p-value
#'  of the above correlation
#' }
#'
#' @examples{
#'  \dontrun{
#'   lesydata = file.path(find.package('LESYMAP'),'extdata')
#'   filenames = Sys.glob(file.path(lesydata, 'lesions', '*.nii.gz'))
#'   behavior = Sys.glob(file.path(lesydata, 'behavior', 'behavior.txt'))
#'   behavior = read.table(behavior,header=FALSE)[,1]
#'   avg = antsAverageImages(filenames)
#'   mask = thresholdImage(avg, 0.1, Inf)
#'   lesmat = imagesToMatrix(filenames,mask)
#'   result = lsm_sccan(lesmat, behavior,
#'        optimizeSparseness=F, sparseness=0.8, mask = mask)
#'  }
#' }
#'
#' @author Dorian Pustina
#'
#' @export
lsm_sccan <- function(lesmat, behavior, mask, showInfo=TRUE,
                      optimizeSparseness = TRUE, validateSparseness=FALSE,
                      pThreshold=0.05,
                      mycoption=1,
                      robust=1,
                      sparseness=0.045,
                      sparseness.behav = -0.99,
                      nvecs=1,
                      cthresh=150,
                      its=20,
                      npermsSCCAN=0,
                      smooth=0.4,
                      maxBased=FALSE,
                      directionalSCCAN=TRUE,
                      ...) {

  sparseness = c( sparseness, sparseness.behav )
  cthresh = c(cthresh,0)

  # scale and center data
  behavior = scale(behavior, scale=T, center=T)
  lesmat = scale(lesmat, scale=T, center=T)
  # prepare data
  inmats=list(lesmat,as.matrix(behavior))
  sccan.masks=c(mask,NA)

  # check if user specified sparseness
  if ('sparseness' %in% names(match.call())) {
    optimizeSparseness = FALSE
    if (!('validateSparseness' %in% names(match.call()))) {
      validateSparseness = TRUE
    }
  }

  if (optimizeSparseness | validateSparseness) {
    sparse.optim = optimize_SCCANsparseness(lesmat = lesmat, behavior = behavior, mask=mask,
                                            cthresh=cthresh, mycoption=mycoption, robust=robust,
                                            nvecs=nvecs, its=its, npermsSCCAN=npermsSCCAN,
                                            smooth=smooth, sparseness.behav=sparseness.behav,
                                            showInfo=showInfo,
                                            maxBased=maxBased,
                                            sparseness=sparseness[1], justValidate=validateSparseness,
                                            directionalSCCAN=directionalSCCAN, ...)

    sparseness = c(sparse.optim$minimum, sparseness.behav)
    CVcorrelation.stat = sparse.optim$CVcorrelation.stat
    r = abs(CVcorrelation.stat)
    n = length(behavior)
    tstat = (r*sqrt(n-2))/(sqrt(1 - r^2))
    CVcorrelation.pval = pt(-abs(tstat), n-2)*2
    CVcorrelation.pval = ifelse(CVcorrelation.pval<1, CVcorrelation.pval, 1) # to fix p > 1

    if (showInfo & !validateSparseness) {
      msg = paste0('\n       Found optimal sparsenes ', round(sparseness[1],3),
                             ' (CV corr=', round(CVcorrelation.stat,3), ' p=', format(CVcorrelation.pval, digits=3), ')')
      printInfo(msg, type='middle')
    }

    if (showInfo & validateSparseness) {
      msg = paste0('\n       Validated sparseness ', round(sparseness[1],3),
                                             ' (CV corr=', round(CVcorrelation.stat,3), ' p=', format(CVcorrelation.pval, digits=3), ')')
      printInfo(msg, type='middle')
    }

    # if poor result, end it here
    if (CVcorrelation.pval > pThreshold) {
      if (showInfo) printInfo('\n       WARNING: Poor cross-validated accuracy, returning NULL result.', type='middle')
      return(list(statistic=rep(0,ncol(lesmat)),
                  pvalue=rep(1,ncol(lesmat)),
                  optimalSparseness = sparse.optim$minimum,
                  CVcorrelation.stat= CVcorrelation.stat,
                  CVcorrelation.pval= CVcorrelation.pval))
    }

  }


  if (showInfo) {
    printInfo(paste('\n       Calling SCCAN with:'))
    printInfo(paste('\n            Components:\t\t', nvecs), type='middle')
    printInfo(paste('\n            Use ranks:\t\t', robust), type='middle')
    printInfo(paste('\n            Sparseness:\t\t', round(sparseness[1], 3)), type='middle')
    printInfo(paste('\n            Cluster threshold:\t', cthresh[1]), type='middle')
    printInfo(paste('\n            Smooth sigma:\t', smooth), type='middle')
    printInfo(paste('\n            Iterations:\t\t', its), type='middle')
    printInfo(paste('\n            maxBased:\t\t', maxBased), type='middle')
    printInfo(paste('\n            directionalSCCAN:\t', directionalSCCAN), type='middle')
    printInfo(paste('\n            optimizeSparseness:\t', optimizeSparseness), type='middle')
    printInfo(paste('\n            validateSparseness:\t', validateSparseness), type='middle')
  }

  sccan = sparseDecom2( inmats,inmask=sccan.masks, mycoption=mycoption,
                           robust=robust, sparseness=sparseness, nvecs=nvecs,
                           cthresh=cthresh,its=its, perms=npermsSCCAN, smooth=smooth,
                        maxBased=maxBased)


  # normalize values to 1 or -1
  statistic = sccan$eig1 / max(abs(sccan$eig1))

  # flip weights if necessary
  if (directionalSCCAN) {
    posbehav = ifelse(sccan$eig2[1,1] < 0, -1, 1)
    poscor = ifelse(sccan$ccasummary$corrs[1] < 0, -1, 1)
    flipval =  posbehav * poscor
    statistic = statistic * flipval
  } else {
    statistic = abs(statistic)
  }

  # shave away weights < 0.1, not needed for maxBased because they are removed in sparseDecom2
  if (!maxBased) statistic[statistic < 0.1 & statistic > -0.1] = 0

  # eleminate small clusters
  # placed on purpose after 0.1 thresholding to remove
  # remaining small leftover clusters
  temp = makeImage(mask,statistic) # put stat in image
  tempclust = labelClusters(abs(temp), minClusterSize = cthresh, minThresh = .Machine$double.eps, maxThresh=Inf)
  temp = temp * thresholdImage(tempclust, .Machine$double.eps, Inf)
  statistic = imageListToMatrix(list(temp), mask)[1,]
  if (showInfo & sum(statistic!=0) == 0) printInfo('\n       WARNING: Post-sccan cluster thresholding removed all voxels.', type='middle')


  output = list(statistic=statistic)

  output$rawWeights.img = makeImage(mask,sccan$eig1)
  output$sccan.eig2 = sccan$eig2
  output$sccan.ccasummary = sccan$ccasummary

  if (optimizeSparseness) {
    output$optimalSparseness = sparse.optim$minimum
    output$CVcorrelation.stat = CVcorrelation.stat
    output$CVcorrelation.pval = CVcorrelation.pval
  }

  return(output)
}
