#' @title Check headers of files stored on disk
#'
#' @description
#' Function to check that all filenames in a vector
#' point to existing files with the same resolution,
#' orientation, size, and origin.
#'
#' @param files character vector of filenames
#' @param showError logical whether to show an error (True) or
#'                    to return a boolean instead. Returned
#'                    values are True=pass,False=Fail
#' @return logical if the test was successful or not
#'
#' @author Dorian Pustina
#'
#' @export
checkFilenameHeaders <- function(files, showError=TRUE) {

  compareto = antsImageHeaderInfo(files[1])
  compnames = names(compareto)
  for (i in 1:length(files)) {
    header = antsImageHeaderInfo(files[i])
    for (h in 1:length(compareto)) {
      if (any(compareto[[h]] != header[[h]])) {
        if (!showError) return(F)
        else stop(paste('File with different', compnames[h],'detected:', files[i]))
      }
    }
  }

  if (!showError) return(TRUE) # if all turned ok
}
