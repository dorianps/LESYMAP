#' @title Compare headers between mask and other images
#'
#' @description
#' Function to check if mask is in the same space as inputs
#'
#' @param lesions.list list of antsImages or character vector
#' of filenames
#' @param mask antsImage of mask to check
#'
#' @return Nothing is returned, function stops with error if
#' mask is not in the same space as images in lesions.list
#'
#' @author Dorian Pustina
#'
#' @export
checkMask <- function(lesions.list, mask) {

  inputtype = checkAntsInput(lesions.list)
  inputtype.mask = checkAntsInput(mask)

  if (inputtype.mask == 'antsFiles' & inputtype == 'antsFiles') { # MASK AND LESIONS ARE FILES
    if (! checkFilenameHeaders(c(mask,lesions.list), showError = F) ) stop('Mask and image(s) are in different space')
#     checkAntsInput(c(mask,lesions.list), checkHeaders = T)
  } else if (inputtype.mask == 'antsImage' & inputtype == 'antsImageList') { # MASK IS IMAGE, LESIONS ARE LIST
    if (!checkImageList(c(mask,lesions.list), showError = F) ) stop('Mask and images are in different space')
#     checkAntsInput(c(mask,lesions.list), checkHeaders = T)
  } else  if (inputtype.mask == 'antsImage' & inputtype == 'antsImage') { # MASK IS IMAGE, LESIONS ARE 4D
    if ( any(antsGetDirection(mask) != antsGetDirection(lesions.list)[1:3,1:3]) ) stop('Mismatch directions of mask and 4D input.')
    if ( any(antsGetOrigin(mask) != antsGetOrigin(lesions.list)[1:3]) ) stop('Mismatch origin of mask and 4D input.')
    if ( any(antsGetSpacing(mask) != antsGetSpacing(lesions.list)[1:3]) ) stop('Mismatch resolution of mask and 4D input.')
    if ( any(dim(mask) != dim(lesions.list)[1:3]) ) stop('Mismatch sizes of mask and 4D input.')
  }

}
