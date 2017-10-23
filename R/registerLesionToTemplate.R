#' registerLesionToTemplate
#'
#' Brings lesion maps in template space by registering the
#' subject\'s anatomical to the template and applying the
#' same transform to the lesion. To improve the registration
#' the anatomical image is bias corrected and denoised. In
#' addition, you can choose to skull-strip the image and run a
#' more careful registration brain-on-brain so that the skull
#' does not impact the registration in any way. Note, for technical
#' reasons the registration is performed counterintuitively by
#' moving the template on the subject, and not the subject on the
#' template. For this reason, to bring the subject in template space
#' we use the inverse transformation. Also note, at the moment ANTsR
#' does not produce an inverse affine transformation explicitly, both
#' forward and inverse affine transforms are identical. You can use ANTs
#' to compute the inverse, or tell ANTsR if you need to invert an affine
#' matrix applying the transformations (see \code{whichtoinvert}).
#'
#' @param subImg antsImage or character filename of the anatomical
#'               image of the subject. Typically this is a T1-weighted
#'               MRI image, on which you drew the lesion map.
#'
#' @param subLesion antsImage or character filename of the lesion map.
#'               Typically you draw this manually or obtain it from
#'               automated lesion segmentation software. You can try
#'               our \href{https://github.com/dorianps/LINDA}{LINDA
#'               toolbox} for an automated alternative. Yet, manual
#'               drawing can be performed
#'               \href{https://youtu.be/ZVmINdWk5R4}{quickly}
#'                and is preferred.
#'
#' @param templateImg antsImage or filename of the anatomical template
#'               image. This image should be with skull included.
#'
#' @param templateBrainMask antsImage or filename of the template brain
#'               mask. This mask is needed for skull-stripped registrations.
#'
#' @param skullStrip logical whether to remove the skull and perform
#'               brain-on-brain registration.
#'
#' @param typeofTransform an \code{antsRegistration} parameter that controls
#'               the quality of registration. The default is \code{SyNCC}, which
#'               probably is the most robust and takes long (1-2 hours maybe).
#'               For faster registration you can try \code{SyN}.
#'
#' @param outprefix character of the prefix where to save the output. If
#'               this is set, most of images and transformations will be
#'               saved at the specified path/prefix. It is passed to
#'               \code{antsRegistration}.
#'
#' @param tstamp format of the timestamp when displaying info messages.
#'
#' @param showInfo logical whether to show info messages or be completely
#'               quiet. If you want also verbose registration messages,
#'               please set \code{verbose=TRUE}.
#'
#' @param ... other arguments to pass to \code{antsRegistration}
#'
#' @export
#' 
#' @return
#' List of objects returned:
#' \itemize{
#'  \item\code{subImg} - subject\'s image in native space (after some preprocessing)
#'  \item\code{subLesion} - subject\'s lesion map in native space
#'  \item\code{subImgTemplate} - subject\'s image in template space
#'  \item\code{subLesionTemplate} - subject\'s lesion in template space
#'  \item\code{registration$inverse_subject2template} - transformation matrices subject to template
#'  \item\code{registration$forward_template2subject} - transformation matrices tempalte to subject
#' }
#'
#' @author Dorian Pustina
#'
#' 
#'

registerLesionToTemplate <- function(subImg, subLesion,
                                     templateImg = NA, templateBrainMask = NA,
                                     skullStrip=T,
                                     typeofTransform = 'SyNCC',
                                     outprefix = '',
                                     tstamp = "%H:%M:%S",
                                     showInfo = T,
                                      ...) {

  
  # load provided MNI template if user did not specify
  if (is.na(templateImg)) {
    if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Template undefined, using provided MNI152_2009c...\n'))
    mnipath = file.path(find.package('LESYMAP'), 'extdata', 'template','other_templates', 'MNI152_2009c')
    templateImg = antsImageRead(file.path(mnipath,'mni_icbm152_t1_tal_nlin_sym_09c.nii.gz'))
    templateBrainMask = antsImageRead(file.path(mnipath,'mni_icbm152_t1_tal_nlin_sym_09c_mask.nii.gz'))
    checkMask(templateImg, templateBrainMask)
  }
  
  
  # load images if filenames are passed
  if ( checkAntsInput(subImg) == 'antsFiles') {
    if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Loading subject\'s anatomical...\n'))
    subImg = antsImageRead(subImg)
  }
  if ( checkAntsInput(subLesion) == 'antsFiles') {
    if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Loading subject\'s lesion file...\n'))
    subLesion = antsImageRead(subLesion)
  }
  if ( !is.na(templateImg) & checkAntsInput(templateImg) == 'antsFiles') {
    if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Loading template anatomical...\n'))
    templateImg = antsImageRead(templateImg)
  }
  if ( !is.na(templateBrainMask) & checkAntsInput(templateBrainMask) == 'antsFiles') {
    if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Loading template brain mask...\n'))
    templateBrainMask = antsImageRead(templateBrainMask)
  }


  # binarize mask inputs to be always 0-1
  subLesion = thresholdImage(subLesion, 0.1, Inf)
  templateBrainMask = thresholdImage(templateBrainMask, 0.1, Inf)

  # make sure lesion and anatomical are in the same space
  if (showInfo)
    cat(paste(format(Sys.time(), tstamp) , 'Assuring lesion and antomical are in the same space...\n'))
  checkMask(subImg, subLesion)


  # truncate outlier intensities, bias correct, and denoise
  if (showInfo)
    cat(paste(format(Sys.time(), tstamp) , 'Running bias correction on anatomical...\n'))
  subImg = abpN4(subImg)
  if (showInfo)
    cat(paste(format(Sys.time(), tstamp) , 'Denoising the anatomical...\n'))
  subImg = subImg %>% iMath('PeronaMalik', 10, 0.4)


  #' here is the rationale
  #' after skull stripping we take the mask and dilate it a bit
  #' then we remove the lesion from it.
  #' At the registration we use the brain-only image coupled with the
  #' mask that is slightly larger. The reason for it is because we want
  #' edge information (brain-air boundary) to improve the registration.
  #' The approach is based on Nick Tustison's advice, see ANTsR issue:
  #' https://github.com/ANTsX/ANTs/issues/483
  #'
  if (skullStrip) {
    if (showInfo)
      cat(paste(format(Sys.time(), tstamp) , 'SKull-stripping subject\'s image...\n'))
    temp = abpBrainExtraction(img = subImg, tem = templateImg, temmask = templateBrainMask,
                              regtype = 'SyNabp')
    subImg = temp$brain
    subBrainMask = temp$bmask
    subBrainMask = subBrainMask %>% iMath('MD', 2)
    subRegMask = subBrainMask - subLesion
    templateImg = templateImg * templateBrainMask
    rm(temp)
  } else {
    subRegMask = subImg * 0 + 1 - subLesion
  }



  # run registration
  # subject needs to be fixed and template moving
  # otherwise can't use mask with existing ANTsR
  if (showInfo) {
    cat(paste(format(Sys.time(), tstamp) , 'Running template registration',
              ifelse(typeofTransform=='SyNCC','(expect >1 hour)',''),'...\n'))
  }
  reg = antsRegistration(fixed = subImg, moving = templateImg, mask = subRegMask,
                   outprefix = outprefix,
                   typeofTransform = typeofTransform)
  # save image if requested
  if (nchar(outprefix) > 0)
    antsImageWrite(reg$warpedfixout, filename = paste0(outprefix, '_subImgTemplate.nii.gz'))

  # bring lesion in template space
  if (showInfo) {
    cat(paste(format(Sys.time(), tstamp) , 'Applying registration to lesion image...\n'))
  }
  subLesionTemplate = antsApplyTransforms(fixed = templateImg, moving = subLesion,
                                      transformlist = reg$invtransforms, whichtoinvert = c(1,0),
                                      interpolator = 'nearestNeighbor')
  # save lesion if requested
  if (nchar(outprefix) > 0)
    antsImageWrite(subLesionTemplate, filename = paste0(outprefix, '_subLesionTemplate.nii.gz'))

  if (showInfo) {
    cat(paste(format(Sys.time(), tstamp) , 'Lesion size native:',
              round( sum(subLesion) * prod(antsGetSpacing(subLesion))/1000, 2),'ml\n'))
    cat(paste(format(Sys.time(), tstamp) , 'Lesion size template:',
              round( sum(subLesion) * prod(antsGetSpacing(subLesionTemplate))/1000, 2),'ml\n'))
  }

  # prepare output list
  output = list()
  output$subImg = subImg
  output$subLesion = subLesion
  output$subImgTemplate = reg$warpedfixout
  output$subLesionTemplate = subLesionTemplate
  output$registration = list()
  output$registration$inverse_subject2template = reg$invtransforms
  output$registration$forward_template2subject = reg$fwdtransforms


  return(output)

}