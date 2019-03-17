#' lesyload_mricron
#'
#' Function to load data from a previous analysis
#' in MRIcron/npm in a ready format for use in lesymap
#'
#' @param valfile mricron filename with extention *.val. The function
#'  will search for images in the same folder where valfile
#'  is located, unless you specify \code{imageFolder}.
#'  If any of the files listed in the .val file are not
#'  found in the folder, an error will be displayed.
#' @param imageFolder (default=NA) folder to look for the image files
#' @param returnFilenames (default=FALSE) By default the function will
#'  load the images in memory to speed up things in lesymap.
#'  This may require too much RAM memory in some cases, and
#'  you may want to use filenames instead, which requires less
#'  memory but is slower in lesymap.
#' @param checkHeaders (default=TRUE) Headers will be checked to make sure
#'  all images have the same dimension/origin/resolution, etc.
#' @param showInfo (default=TRUE) show information upon successful load
#'
#' @return List with the following information
#'      lesions - list of antsImages or vector of filenames
#'      behavior - vector of behavioral scores
#'
#' @author Dorian Pustina
#'
#' @export
lesyload_mricron <- function(valfile, imageFolder=NA, returnFilenames=F, checkHeaders=T, showInfo=T) {

  valdata = read.table(valfile, header = T, sep = '\t', dec = '.', strip.white = T, stringsAsFactors = F)

  if (is.na(imageFolder)) imageFolder = dirname(valfile)
  valdata$ImageName = file.path(imageFolder, valdata$ImageName)

  # check files exist
  existindx = file.exists(valdata$ImageName)
  if (sum(!existindx) == nrow(valdata)) {
    stop('\n  None of the image files are found.\n  You may want to specify their location with imageFolder.\n')
  } else if (any(!existindx)) {
    stop(paste0('  Some of the files are not found:\n', paste(basename(valdata$ImageName[!existindx]), collapse='\n')))
  }

  if (showInfo) printInfo('All files found.')

  # all files are there, check their headers
  if (checkHeaders) {
    if (showInfo) printInfo('Checking for unusual images...', type='head')
    checkFilenameHeaders(files = valdata$ImageName, showError = F)
    if (showInfo) printInfo('all headers identical.', type='tail')
  }

  # prepare lesions output
  if (returnFilenames) {
    if (showInfo) printInfo('Filenames will be returned.')
    lesions = valdata$ImageName
  } else {
    if (showInfo) printInfo('Image list will be returned, please wait...', type='head')
    lesions = imageFileNames2ImageList(valdata$ImageName)
    if (showInfo) printInfo('all loaded.', type='tail')
  }

  # prepare behavior output
  behavior = valdata[,2]

  if (showInfo) printInfo('You can now run lesymap with these data.')

  invisible(gc())
  return(list(lesions=lesions,
              behavior=behavior))

}
