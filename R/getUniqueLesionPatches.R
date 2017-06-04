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
#'  \itemize{
#'  \item\code{patchimg} - antsImage with every voxel assigned a patch number
#'  \item\code{patchimg.samples} - antsImage mask of one representative voxel
#'  for each patch. Can be used to extract the patchmatrix.
#'  \item\code{patchimg.size} - antsImage with the patch size at every voxel
#'  \item\code{patchimg.mask} - antsImage of full mask with all voxels. Can
#'  be used to put back results in combination with \code{patchindx}.
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
#' \dontrun{
#' onefile = system.file(file.path('extdata','lesions'), 'Subject_001.nii.gz', package='LESYMAP')
#' niftifolder = dirname(onefile)
#'
#' filenames = Sys.glob( file.path(niftifolder, 'Subject*.nii.gz'))
#' patches = getUniqueLesionPatches(filenames) # slower
#'
#' lesions = imageFileNames2ImageList(filenames)
#' patches = getUniqueLesionPatches(lesions) # faster
#' }
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
