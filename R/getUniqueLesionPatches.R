#' Unique Lesion Patches
#'
#' Compute uniqe patches of voxels with the same
#' pattern of lesions in all subjects. Useful to understand
#' the number of patterns that will be analyzed in a lesion dataset.
#' A patch is a group of voxels, not necessarily close to each
#' other, which have the same identical lesion pattern.
#'
#' @param lesions.list list of antsImages (faster) or filenames (slower)
#' @param mask (default=NA) a mask image to restrict the search for patches. Will be
#' automatically calculated if not provided. Normally the mask restricts
#' the search only to voxels lesioned in >10\% of subejcts. To set this proportion
#' use \code{thresholdPercent}.
#' @param returnPatchMatrix (default=FALSE) logical, should the matrix of patches
#' be returned. This is used in \code{\link{lesymap}} to run the analyses.
#' @param thresholdPercent (default=0.1) voxels with lesions in less than
#' this proportion of subjects will not be considered. I.e., 0.1 = 10\%.
#' @param binaryCheck (default=FALSE) set this to TRUE to verify that maps are binary.
#' @param showInfo (default=TRUE) logical indicating whether to display information.
#'
#' @return
#'  List of objects named as follows:
#'  \itemize{
#'  \item\code{patchimg} - antsImage with every voxel assigned a patch number
#'  \item\code{patchimg.samples} - antsImage mask of one representative voxel
#'  for each patch. Can be used to extract the patchmatrix.
#'  \item\code{patchimg.size} - antsImage with the patch size at every voxel
#'  \item\code{patchimg.mask} - antsImage of the mask used to extract patches. Can
#'  be used to put back results when combined with \code{patchindx}.
#'  \item\code{patchindx} - vector of patch membership for each voxel. Can
#'  be used to put back results in an image.
#'  \item\code{npatches} - number of unique patches in the image
#'  \item\code{nvoxels} - total number of lesioned voxels in \code{patchimg.mask}
#'  \item\code{patchvoxels} - vector of voxel count for each patch
#'  \item\code{patchvolumes} - vector of volume size for each patch
#'  \item\code{patchmatrix}  - matrix of patches. This is used in lesymap to
#'  save time when running repetitive analyses.
#'  }
#'
#'
#' @examples
#' lesydata = file.path(find.package('LESYMAP'),'extdata')
#'
#' filenames = Sys.glob(file.path(lesydata, 'lesions', '*.nii.gz'))
#' patchinfo = getUniqueLesionPatches(filenames[1:10]) # slower
#'
#' lesions = imageFileNames2ImageList(filenames[1:10])
#' patchinfo = getUniqueLesionPatches(lesions) # faster
#'
#' @author Dorian Pustina
#'
#' @export
getUniqueLesionPatches <- function(lesions.list, mask=NA, returnPatchMatrix=F,
                                       thresholdPercent=0.1, binaryCheck=F, showInfo=T) {

  # cannot set both mask and thresholdPercent
  if (!is.na(mask) & thresholdPercent!=0.1) stop('Cannot set both mask and thresholdPercent. Choose one.')
  if (thresholdPercent < 0 | thresholdPercent > 1) stop('thresholdPercent must be between 0 and 1.')

  inputtype = checkAntsInput(lesions.list, checkHeaders = F)

  # compute mask from average, if not defined by user
  if (is.na(mask)) {
    if (inputtype == 'antsImage') avgles = getAverageOfTimeSeries(lesions.list)
    else avgles = antsAverageImages(lesions.list)

    if (max(avgles) > 1) stop('Average lesion map exceeds 1. This is not possible.')
    # compute mask, no too few or too many subjects lesioned
    mask = thresholdImage(avgles, thresholdPercent, 1 - thresholdPercent)
    if (thresholdPercent == 0) mask[avgles==0] = 0 # remove 0 voxels from mask
  }

  # put images to matrix
  if (inputtype == 'antsImageList') { # list of antsImages
    checkImageList(lesions.list, binaryCheck = binaryCheck)
    lesmat = imageListToMatrix(lesions.list, mask)
  } else if (inputtype == 'antsFiles') { # vector of filenames
    lesmat = imagesToMatrix(lesions.list, mask)
  } else if (inputtype == 'antsImage') {
    lesmat = timeseries2matrix(lesions.list, mask)
  }

  # engine for of unique patch computation
  add = 1
  summed = rep(0, ncol(lesmat))
  for (i in 1:nrow(lesmat)) {
    summed = summed + (lesmat[i, ]*add)
    summed = match(summed, unique(summed))
    add = max(summed)+1
  }



  # compute statistics
  numpatches = length(unique(summed))
  numvoxels = ncol(lesmat)

  # put patch numbers in the image
  patchimg = mask*0
  patchimg[mask==1] = summed

  patchstat = labelStats(patchimg,patchimg)
  if (any(patchstat$LabelValue %in% 0)) patchstat = patchstat[ !(patchstat$LabelValue %in% 0) ,]

  # compute image with patch size values
  patchsize = patchimg*0
  patchsize[mask==1] = patchstat$Count[summed]

  # compute patchsamples
  patindx = as.numeric(!duplicated(summed))
  patchsamples = patchimg*0
  patchsamples[mask==1] = patindx

  # compute patch matrix if requested
  if (returnPatchMatrix) patchmatrix = lesmat[ , patindx==1]

  if (showInfo) {
    printInfo(paste(
      numpatches, 'unique patches,',
      numvoxels, 'voxels -',
      round(numvoxels/numpatches,1), 'times more voxels'
    ), type='tail')
  }

  output = list(
    patchimg=patchimg,
    patchimg.samples=patchsamples,
    patchimg.size=patchsize,
    patchimg.mask=mask,
    patchindx=summed,
    npatches=numpatches,
    nvoxels=numvoxels,
    patchvoxels=patchstat$Count,
    patchvolumes=patchstat$Volume)

  if (returnPatchMatrix) output$patchmatrix = patchmatrix

  invisible(gc())
  return(output)
}
