#' getLesionSize
#'
#' Compute lesion sizes from a list of antsImages.
#'
#' @param lesions.list List of antsImages or vector of filenames.
#'                     It is assumed that images are binary (0/1).
#' @param showInfo logical show or not informations/warnings
#'
#' @return vector of lesion sizes in mm3
#'
#' @author Dorian Pustina
#'
#' @export
getLesionSize <- function(lesions.list, showInfo=TRUE) {

  inputtype = checkAntsInput(lesions.list, checkHeaders = F)
  if ( !(inputtype %in% c('antsImageList', 'antsFiles', 'antsImage')) ) {
    stop('Unrecognized input: lesions.list can be a 4D image, a list of images, or a vector of filenames')
  }

  if (inputtype == 'antsImage') {
    lesions.list = list(lesions.list)
    inputtype = 'antsImageList'
  }

  if (inputtype == 'antsFiles' & showInfo)
    printInfo('      Slow computation from filenames. Preload images with imageFileNames2ImageList() for faster processing.', type='middle')

  lessize = rep(NA, length(lesions.list))

  for (i in 1:length(lesions.list)) {

    if (inputtype == 'antsFiles') {
      thisimg = antsImageRead(lesions.list[i])
      invisible(gc()) # release memory to help R in linux
    } else {
      thisimg = lesions.list[[i]]
    }

    voxsize = prod(antsGetSpacing(thisimg))
    # strangely faster with as.array
    # https://github.com/ANTsX/ANTsR/issues/233
    # may need to change after ANTsR is fixed
    lessize[i] = sum(as.array(thisimg)) * voxsize

  }

  return(lessize)
}
