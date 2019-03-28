#' @title Check headers of list of antsImages
#'
#' @description
#' Function to check that all antsImages in a list have
#' the same orientation, origin, and resolution.
#' The function stops with an error if one of the images
#' has unusual headers. This behavior can be overcome
#' by setting showError=F, and using the returned status
#' (True=pass, False=fail) to make decisions outside this
#' function.
#'
#' @param imgList list of antsImages
#' @param showError boolean indicating whether to show the
#' exact error and interrupt the function (TRUE, default),
#' or don't show the error and return the check
#' status (FALSE). The returned values when
#' showError=F are T=passed or F=Failed.
#' @param binaryCheck boolean, check if images are binary (0/1 values).
#' Useful when checking masks or lesions.
#' This check slows the output of the function.
#'
#' @return True if list has images with same headers,
#'          otherwise False.
#'
#' @examples
#' \dontrun{
#' files = Sys.glob('/data/jag/nifti/*.nii.gz')
#' myimagelist = imageFileNames2ImageList(files)
#' checkImageList(myimagelist) # no value returned
#' checkImageList(lesions, showError=F) # True returned
#' myimagelist[[4]] = cropIndices(myimagelist[[4]], c(1,1,1), c(20,20,20))
#' checkImageList(myimagelist) # error on image 4
#' }
#'
#' @author Dorian Pustina
#'
#'@export
checkImageList <- function(imgList, showError=T, binaryCheck=F) {

  # check this is a list
  if (class(imgList) != 'list')  {
    if (!showError) return(F)
    else stop('Input is not a list of antsImages.')
  }

  # check all images are antsImages
  for (i in 1:length(imgList)) if (class(imgList[[i]]) != 'antsImage') {
    if (!showError) return(F)
    else stop(paste('Image',i,'is not antsImage'))
  }

  # check all images have same:
  # resolution
  compareto = antsGetSpacing(imgList[[1]])
  for (i in 1:length(imgList)) if (any(compareto != antsGetSpacing(imgList[[i]]))) {
    if (!showError) return(F)
    else stop(paste('Image with different dimension detected:', i))
  }
  # orientation
  compareto = antsGetDirection(imgList[[1]])
  for (i in 1:length(imgList)) if (any(compareto != antsGetDirection(imgList[[i]]))) {
    if (!showError) return(F)
    else stop(paste('Image with different orientation detected:', i))
  }
  # origin
  compareto = antsGetOrigin(imgList[[1]])
  for (i in 1:length(imgList)) if (any(compareto != antsGetOrigin(imgList[[i]]))) {
    if (!showError) return(F)
    else stop(paste('Image with different origin detected:', i))
  }
  # dimensions
  compareto = dim(imgList[[1]])
  for (i in 1:length(imgList)) if (any(compareto != dim(imgList[[i]]))) {
    if (!showError) return(F)
    else stop(paste('Image with different dimension detected:', i))
  }

  # check image is binary
  if (binaryCheck) {
    for (i in 1:length(imgList)) if (any(!unique(c(as.numeric(imgList[[i]]))) %in% c(0,1) )) {
      if (!showError) return(F)
      else stop(paste('Non binary image', i))
    }
  }

  if (!showError) return(T)
}
