#' @title Prediction of new cases from lesymap output
#'
#' @description
#' Uses an existing lesyamp object output from your
#' analysis to predict new cases.
#'
#' @param lsm object of class lesymap from previous analysis
#'
#' @param lesions.list list of antsImages, or a vector of
#' filenames, or a single antsImage with 4 dimensions.
#'
#' @param binaryCheck logical (default=FALSE), make sure the
#' lesion matrix is 0/1. This will help if lesion maps are drawn
#' in MRIcron or other software which label lesioned voxel
#' with value 255.
#'
#' @param showInfo logical (default=TRUE), display time-stamped info messages
#'
#' @param ... other arguments for flexible calling from other functions.
#'
#' @return
#' Vector of predicted values:
#' \itemize{
#'  \item\code{behavior.scaled} - scaled values as predicted by the model
#'  \item\code{behavior.raw} - descaled raw values
#' }
#'
#' @examples{
#'  \dontrun{
#'   lesydata = file.path(find.package('LESYMAP'),'extdata')
#'   filenames = Sys.glob(file.path(lesydata, 'lesions', 'Subject*.nii.gz'))
#'   behavior = Sys.glob(file.path(lesydata, 'behavior', 'behavior.txt'))
#'   lesions = imageFileNames2ImageList(filenames)
#'   behav = read.table(behavior)$V1 * 1000
#'
#'   train = 1:100
#'   test = 101:131
#'
#'   lsm = lesymap(lesions[train], behav[train], method='sccan',
#'   sparseness=0.2, validateSparseness=F)
#'   predbehav = lesymap.predict(lsm, lesions[test])
#'  }
#' }
#'
#' @author Dorian Pustina
#'
#' @export

lesymap.predict <- function(lsm, lesions.list,
                            binaryCheck = TRUE,
                            showInfo = TRUE,
                            ...) {

  # make sure this is a lesymap object
  if (! 'lesymap' %in% class(lsm)) stop('lsm must be an object of class lesymap')

  # make sure we can predict this method
  if (! lsm$callinfo$method %in% c('sccan') ) {
    stop(paste('Analysis with <', lsm$callinfo$method,
               '> cannot be predicted, only < sccan > is currently enabled for prediction.'))
  }


  ### check input type
  inputtype = checkAntsInput(lesions.list, checkHeaders = TRUE)

  # make sure input is a list or filenames
  if ( !(inputtype %in% c('antsImageList', 'antsFiles', 'antsImage')) ) {
    stop('Unrecognized input: lesions.list should be a 4D antsImage, a list of antsImages, or a vector of filenames')
  }

  ### start dealing with lesions.list
  ############## 4D case
  # single filename, check it's 4D, load it if filename, unpack it to list
  if (inputtype == 'antsFiles' & length(lesions.list) == 1) {
    temp = antsImageHeaderInfo(lesions.list[1])
    if (temp$nDimensions != 4) stop('File is not a 4D image. You must point to a 4D file when a single filename is defined.')
    if (showInfo) printInfo('Loading 4D image...', type='head')
    lesions.list = antsImageRead(lesions.list[1])
    if (showInfo) printInfo(paste( dim(lesions.list)[4], 'images present.'), type='tail')
    if (showInfo) printInfo('Converting 4D image into 3D image list...')
    lesions.list = splitNDImageToList(lesions.list)
    invisible(gc()) # free some memory after conversion
    inputtype = 'antsImageList'
  } else if (inputtype == 'antsImage') {
    if (lesions.list@dimension != 4) stop('Input is a single image but is not 4D.')
    if (showInfo) printInfo('Single 4D image passed, converting to 3D image list...')
    lesions.list = splitNDImageToList(lesions.list)
    invisible(gc()) # free some memory after conversion
    inputtype = 'antsImageList'
  }


  ##############
  # Few tests on images coming as filenames
  # check proper binarization and 255 values, maybe preload
  if (inputtype == 'antsFiles') {
    if (showInfo) printInfo('Filenames as input, checking lesion values on 1st image...')

    temp = antsImageRead(lesions.list[1])
    voxvals = unique(c(as.numeric(temp)))

    if (length(voxvals) != 2) stop('Non binary image detected. Lesions should have only two values (0/1).')

    if (any(! voxvals %in% c(0,1) )) {
      if (showInfo) printInfo('Detected unusual lesion values, loading files into memory to fix...')
      lesions.list = imageFileNames2ImageList(lesions.list)
      inputtype = 'antsImageList'
    }
  }



  ##############
  # image list might be from MRIcron, convert to binary
  # if input='antsFiles', it needs a binary check later on lesmat
  if (inputtype == 'antsImageList') {
    rebinarize = FALSE
    if (max(lesions.list[[1]]) > 1) rebinarize = TRUE # just check 1st, for is too long
    # for (i in 1:length(lesions.list)) {
    #   if (max(as.array(lesions.list[[i]])) > 1) {
    #     rebinarize = TRUE
    #     break
    #   }
    # }

    # perform binarization if needed
    if (rebinarize) {
      if (showInfo) printInfo('Detected lesion value above 1. Rebinarizing 0/1...')
      for (i in 1:length(lesions.list)) lesions.list[[i]] = thresholdImage(lesions.list[[i]], 0.1, Inf)
      binaryCheck = FALSE # no need to check binarization anymore
    }
  }


  #########
  # check antsImageList is binary
  # for antsFiles, we checked only 1st, and
  # will check lesmat later
  if (binaryCheck & inputtype == 'antsImageList') {
    if (showInfo) printInfo('Verifying that lesions are binary 0/1...')
    checkImageList(lesions.list, binaryCheck = TRUE)
  }

  ###
  # at some point we must check headers of lesions.list
  # are same as mask in lsm object
  ###


  #######
  # we are ready to create the voxel matrix
  voxmask = lsm$mask.img

  if (showInfo) printInfo('Computing lesion matrix... ', type='head')
  if (inputtype == 'antsImageList') { # list of antsImages
    lesmat = imageListToMatrix(lesions.list, voxmask)
  } else if (inputtype == 'antsFiles') { # vector of filenames
    lesmat = imagesToMatrix(lesions.list, voxmask)
  }

  if (showInfo) printInfo(paste(dim(lesmat), collapse='x'), type='tail')
  ### we got lesmat now


  #######
  # load weights
  weights = imageListToMatrix(list(lsm$rawWeights.img), voxmask)
  ### got weghts now

  #######
  # get eig2
  eig2 = lsm$sccan.eig2
  #######


  #########
  # predict behavior
  lesmat = scale(lesmat, scale=lsm$sccan.lesmat.scaleval, center=lsm$sccan.lesmat.centerval)
  behavior.scaled = lesmat %*% t(weights) %*% eig2
  #########

  ########
  # unscale behavior
  behavior.raw = behavior.scaled * lsm$sccan.behavior.scaleval + lsm$sccan.behavior.centerval
  behavior.raw = predict(lsm$sccan.predictlm, newdata = data.frame(predbehav.raw = behavior.raw))
  ########

  output = list()
  output$behavior.scaled = behavior.scaled
  output$behavior.raw = behavior.raw
  return(output)
}
