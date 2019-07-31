#' @title Compute regional lesion load
#'
#' @description
#' Computes lesion loads from a series of images. A parcellation
#' image (or simple mask) is required to define the regions from
#' which to compute the lesion load.
#'
#' @param lesions.list list of antsImages or filenames.
#'  Must be binary (0 and 1 values).
#' @param parcellation ansImage or filename of the parcellated volumes. A
#' parcellation is an image brain regions showned as with integer values
#' (i.e. ,1,2,3,...).
#' @param label (default=NA) you can ask to get output for a specific label
#' in the parcellation volume (i.e., label=122).
#' @param mask (default=NA) if this mask is specified (antsImage or filename)
#' lesioned voxels outside the mask are ignored. This is not a
#' good choice, but in case you need it its there.
#' @param binaryCheck (default=FALSE) check whether lesion maps are binary (0/1).
#' Will output an error if lesion files are not binary.
#' @param keepAllLabels (default=FALSE) by default labels are removed if affected
#' in just few subjects. Setting this to TRUE will keep all labels.
#' @param minSubjectPerLabel minimum number of subjects a parcel
#' must be lesioned to keep and return it.
#'
#' @return
#' \itemize{
#'  \item{output}{Matrix of lesion loads between 0 and 1. 1 means 100\% lesioned.
#'  Each column is a single parcel and each row a single subject. Parcel
#'  numbers are placed as column names.}
#'  }
#'
#' @examples
#' lesydata = file.path(find.package('LESYMAP'),'extdata')
#' filenames = Sys.glob(file.path(lesydata, 'lesions', '*.nii.gz'))
#' lesions = imageFileNames2ImageList(filenames[1:10])
#' parcellation = antsImageRead(
#' file.path(lesydata,'template', 'Parcellation_403areas.nii.gz'))
#' lesload = getLesionLoad(lesions, parcellation)
#'
#' @author Dorian Pustina
#'
#' @export
getLesionLoad <- function(lesions.list, parcellation, label=NA,
                          mask=NA, binaryCheck=F, keepAllLabels=F,
                          minSubjectPerLabel = '10%') {

  inputtype = checkAntsInput(lesions.list)

  ###### RUN CHECKS ######

  # make sure input is a list or filenames
  if ( !(inputtype %in% c('antsImageList', 'antsFiles')) ) {
    stop('VLSM input must be a list of antsImages or a vector of filenames')
  }

  # check parcellation is antsImage
  if ( !(checkAntsInput(parcellation) %in% c('antsImage', 'antsFiles')) ) {
    stop('Parcellation argument is not an antsImage or a filename.')
  }

  # check parcellation is image if lesions are antsImageList
  if ( (inputtype == 'antsImageList') & (checkAntsInput(parcellation) != 'antsImage') ) {
    stop('Lesions and parcellations should be both filenames or antsImages (1)')
  } else if ((inputtype == 'antsFiles')  & (checkAntsInput(parcellation) != 'antsFiles')) {
    stop('Lesions and parcellations should be both filenames or antsImages (3)')
  }

  # check parcellation has same headers as lesions
  if (inputtype == 'antsImageList') {
    if (! checkImageList(c(parcellation,lesions.list), binaryCheck = F, showError = F))
      stop('Parcellation image does not have the same headers as lesion images')
  } else if (inputtype == 'antsFiles') {
    if (! checkFilenameHeaders(c(parcellation,lesions.list), showError = F))
      stop('Parcellation image does not have the same headers as lesion images')
  }

  # check images are all binary
  if (binaryCheck & inputtype == 'antsImageList') checkImageList(lesions.list, binaryCheck = T)

  # check the mask
  if (!is.na(c(mask))) {
    if (checkAntsInput(mask) == 'antsFiles') {
      if (inputtype == 'antsFiles') if (!checkFilenameHeaders(c(mask,lesions.list), showError = F))
        stop('Mask have different headers from lesions')
      mask=antsImageRead(mask)
    }
    if (inputtype == 'antsImageList') if (!checkImageList(list(mask), showError = F))
      stop('Mask have different headers from lesions in list.')
  }

  # compute thresholdPercent to remove labels with too few subjects
  if (!is.numeric(minSubjectPerLabel) & is.character(minSubjectPerLabel)) { # input is percentage
    thresholdPercent = as.numeric(gsub('%','', minSubjectPerLabel)) / 100
  } else if (is.numeric(minSubjectPerLabel)) { # user defined exact subject number
    thresholdPercent = minSubjectPerLabel / nrow(lesload)
  }


  ###### START THE WORK ######

  # load parcellation image eventually
  if (inputtype == 'antsFiles') parcellation = antsImageRead(parcellation)

  # keep only required labels eventually
  if (all(!is.na(label))) {
    labindx = as.array(parcellation) %in% label
    parcellation[!labindx] = 0
  }

  # compute lesion load matrix
  temp = labelStats(parcellation,parcellation)
  lesload = matrix(NA, nrow=length(lesions.list), ncol=nrow(temp))
  colnames(lesload) = temp$LabelValue
  for (i in 1:length(lesions.list)) {
    # load or select image from list
    if (inputtype == 'antsFiles') {
      invisible(gc()) # need to release memory in linux
      img = antsImageRead(lesions.list[i])
      if (binaryCheck & !checkImageList(list(img), binaryCheck=binaryCheck, showError=F))
        stop(paste('File is not binary:', lesions.list[i]))
    } else {
      img = lesions.list[[i]]
    }

    # multiply lesion with mask eventually
    if (class(mask) == 'antsImage') img = img*mask

    # put labelStats into matrix
    lesload[i,] = labelStats(img, parcellation)$Mean
  }


  # remove labels affected in too few subjects
  if (!keepAllLabels) {
    lesload = lesload[ , !colnames(lesload) %in% 0] # remove label 0
    lesload = as.matrix(lesload)
    if (thresholdPercent > 0) threshindx = (colSums(lesload>0)/nrow(lesload)) >= thresholdPercent
    else threshindx = (colSums(lesload>0)/nrow(lesload)) > 0
    lesload = lesload[ ,threshindx] # remove unaffacted labels
    lesload = as.matrix(lesload)
  }


  return(lesload)
}
