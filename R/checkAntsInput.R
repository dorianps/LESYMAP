#' @title Determine type of variable passed by user
#'
#' @description
#' Function to check a variable whether is composed
#' of an antsImage, list of antsImages, or simply filenames.
#' If none of the above, an error is returned.
#'
#' @param input the variable to be checked
#' @param checkHeaders make sure all images have
#'                        the same headers
#'
#' @return Type of variable (antsImage, antsImageList, antsFiles)
#'          or error if variable cannot be established.
#'
#' @examples
#' \dontrun{
#'         files = Sys.glob('/data/jag/nifti/*.nii.gz')
#'          myimagelist = imageFileNames2ImageList(files)
#'          checkAntsInput(myimagelist) # returns 'antsImageList'
#'          checkAntsInput(antsFiles) # returns 'antsFiles'
#'          checkAntsInput(myimagelist[[1]]) # returns 'antsImage'
#'          }
#' @author Dorian Pustina
#'
#' @export
checkAntsInput <- function(input, checkHeaders=FALSE) {

  if (class(input) == 'antsImage') return('antsImage')
  else if (class(input) == 'list') { # list of antsImages
    # verify all elements are antsImages
    for (i in 1:length(input)) if (class(input[[i]]) != 'antsImage') {
      stop('Elements other than antsImage detected in input list.')
    }
    # list has only antsImages, final check and return
    if (checkHeaders) checkImageList(input)
    return('antsImageList')
  } else if (class(input) == 'character') { # vector of filenames
    if (any(!file.exists(input))) stop('Character input contains inexistent files')
    if (checkHeaders) checkFilenameHeaders(input)
    return('antsFiles')
  } else if (class(input) == 'numeric') {
    return('vector')
  } else {
    stop('Cannot establish input.')
  }

}
