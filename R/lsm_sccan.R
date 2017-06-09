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
#' @param rawStat logical (default=FALSE) whether to skip
#' the normalization of values and thresholding at 0.1.
#' If TRUE, the raw voxel weights will be returned.
#' @param optimizeSparseness logical (default=TRUE) whether to
#' run the sparseness optimization routine. If false, the defau
#' sparseness value will be used. If sparseness is manually set
#' this flag decides if the manual sparseness will be checked
#' with cross validations.
#' @param pThreshold (default=0.05) If cross validated
#' correlations show significance below this value
#' the results are considered null and an empty
#' map is returned.
#' @param showInfo logical (default-TRUE) display messages
#' @param tstamp timestamp format used in LESYMAP
#' @param sparseness (default=1) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param mycoption (default=1) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param robust (ddefault=1) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param nvecs (default=1) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param cthresh (default=150) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param its (default=20) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param smooth (default=0.4) SCCAN parameter, see \code{\link[ANTsR]{sparseDecom2}}
#' @param npermsSCCAN (default=0) SCCAN permutations, see \code{\link[ANTsR]{sparseDecom2}}
#' @param ... other arguments received from \code{\link{lesymap}}.
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{statistic} - vector of statistical values
#' \item\code{pvalue} - vector of pvalues
#' \item\code{optimalSparseness} - (optional) optimal value found for sparseness
#' \item\code{CVcorrelation.stat} - (optional) Correlation between
#' true and predicted score with k-fold validation using
#' the optimal sparseness value
#' \item\code{CVcorrelation.pval} - (optional) p-value of the above correlation
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
lsm_sccan <- function(lesmat, behavior, mask, rawStat=F, showInfo=T, optimizeSparseness = T,
                      tstamp = "%H:%M:%S", pThreshold=0.05,
                      mycoption=1,
                      robust=1,
                      sparseness=0.045,
                      nvecs=1,
                      cthresh=150,
                      its=20,
                      npermsSCCAN=0,
                      smooth=0.4,
                      ...) {

  sparseness.behav = -0.99
  sparseness = c( sparseness, sparseness.behav )
  cthresh = c(cthresh,0)

  # scale and center data
  behavior = scale(behavior, scale=T, center=T)
  lesmat = scale(lesmat, scale=T, center=T)
  # prepare data
  inmats=list(lesmat,as.matrix(behavior))
  sccan.masks=c(mask,NA)

  # check if user specified sparseness
  justValidate = F
  if ('sparseness' %in% names(match.call())) justValidate = T

  if (optimizeSparseness) {
    if (showInfo & !justValidate) cat(paste('\n       Searching for optimal sparseness:'))
    if (showInfo & justValidate) cat(paste('\n       Validating sparseness:'))

    sparse.optim = optimize_SCCANsparseness(lesmat = lesmat, behavior = behavior, mask=mask,
                                            cthresh=cthresh, mycoption=mycoption, robust=robust,
                                            nvecs=nvecs, its=its, npermsSCCAN=npermsSCCAN,
                                            smooth=smooth, sparseness.behav=sparseness.behav,
                                            showInfo=showInfo, tstamp=tstamp,
                                            sparseness=sparseness[1], justValidate=justValidate, ...)

    sparseness = c(sparse.optim$minimum, sparseness.behav)
    CVcorrelation.stat = sparse.optim$CVcorrelation.stat
    r = abs(CVcorrelation.stat)
    n = length(behavior)
    tstat = (r*sqrt(n-2))/(sqrt(1 - r^2))
    CVcorrelation.pval = pt(-abs(tstat), n-2)*2
    CVcorrelation.pval = ifelse(CVcorrelation.pval<1, CVcorrelation.pval, 1) # to fix p > 1

    if (showInfo & !justValidate) cat(paste0('\n       Found optimal sparsenes ', round(sparseness[1],3),
                             ' (CV corr=', round(CVcorrelation.stat,3), ' p=', format(CVcorrelation.pval, digits=3), ')'))

    if (showInfo & justValidate) cat(paste0('\n       Validated sparseness ', round(sparseness[1],3),
                                             ' (CV corr=', round(CVcorrelation.stat,3), ' p=', format(CVcorrelation.pval, digits=3), ')'))

    # if poor result, end it here
    if (CVcorrelation.pval > pThreshold) {
      if (showInfo) warning('\n       Poor cross-validated accuracy, returning NULL result.')
      return(list(statistic=rep(0,ncol(lesmat)),
                  pvalue=rep(1,ncol(lesmat)),
                  optimalSparseness = sparse.optim$minimum,
                  CVcorrelation.stat= CVcorrelation.stat,
                  CVcorrelation.pval= CVcorrelation.pval))
    }

  }


  if (showInfo) {
    cat(paste('\n       Calling SCCAN with:'))
    cat(paste('\n            Components =', nvecs))
    cat(paste('\n            Use ranks =', robust))
    cat(paste('\n            Sparseness =', round(sparseness[1], 3)))
    cat(paste('\n            Cluster threshold =', cthresh[1]))
    cat(paste('\n            Smooth sigma =', smooth))
    cat(paste('\n            Iterations =', its))
    cat(paste('\n            Permutations =', npermsSCCAN))
  }

  sccan = sparseDecom2( inmats,inmask=sccan.masks, mycoption=mycoption,
                           robust=robust, sparseness=sparseness, nvecs=nvecs,
                           cthresh=cthresh,its=its, perms=npermsSCCAN, smooth=smooth )

  if (!rawStat) {
    statistic = sccan$eig1 / max(abs(sccan$eig1))
    statistic = abs(statistic)
    statistic[statistic < 0.1] = 0
    # eleminate small clusters
    temp = makeImage(mask,statistic)
    tempclust = labelClusters(temp, minClusterSize = cthresh, minThresh = .Machine$double.eps, Inf)
    temp = temp * thresholdImage(tempclust, .Machine$double.eps, Inf)
    statistic = imageListToMatrix(list(temp), mask)[1,]
  } else {
    statistic = sccan$eig1
  }
  pvalue = statistic*0

  output = list(statistic=statistic, pvalue=pvalue)

  if (optimizeSparseness) {
    output$optimalSparseness = sparse.optim$minimum
    output$CVcorrelation.stat = CVcorrelation.stat
    output$CVcorrelation.pval = CVcorrelation.pval
  }

  return(output)
}
