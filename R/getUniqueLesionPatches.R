#' Unique Lesion Patches
#'
#' Compute uniqe patches of voxels with the same
#' pattern of lesions in all subjects. Useful to understand
#' the amount of unique information in a lesion dataset.
#' I.e., your dataset can have 1000 subjects but still have
#' little spatial discrimination because all lesions are
#' similar
#'
#' @param lesions.list list of antsImages (faster) or filenames (slower)
#' @param mask a mask image to restrict the search for patches. Will be
#' automatically calculated if not provided (voxels > 0 in at
#' least one subject)
#' @param returnPatchMatrix logical, should the matrix of patches
#' be returned
#' @param thresholdPercent voxels with lesions in less than X percent of
#'  subjects will not be considered (default 10\%)
#' @param binaryCheck logical, if True images will be verified that are binary
#' @param showInfo logical indicating whether to display information (default=T)
#'
#' @return
#'  List of objects named as follows:
#'    patchimg - antsImage with the patch number each voxels belongs to
#'    patchimg.samples - antsImage mask with a single voxel per patch
#'    patchimg.size - antsImage with the patch size at each voxel
#'    patchimg.mask - the mask within which the function looked for patches
#'    npatches - number of unique patches in the image
#'    nvoxels - total number of lesioned voxels, computed as all voxels
#'       lesioned in at least one subject, within the specified mask
#'    patchvoxels - vector of voxel count for each patch
#'    patchvolumes - vector of volume size for each patch
#'    patchmatrix  - matrix of patches
#'
#'
#'  @examples
#'  \dontrun{
#'  files = Sys.glob('/data/jag/VLSM/*.nii.gz')
#'  lesions = imageFileNames2ImageList(files)
#'
#'  patches = getUniqueLesionPatches(files) # slower
#'  patches = getUniqueLesionPatches(lesions) # faster
#'  avgles = antsAverageImages(lesions)
#'  mask = thresholdImage(avgles, 0.1, Inf)
#'  patches = getUniqueLesionPatches(lesions, mask=mask) # even faster
#'  }
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
    # compute mask
    mask = thresholdImage(avgles, thresholdPercent, Inf)
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
    cat(paste(
      numpatches, 'unique patches,',
      numvoxels, 'voxels -',
      round(numvoxels/numpatches,1), 'times more voxels\n'
    ))
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
