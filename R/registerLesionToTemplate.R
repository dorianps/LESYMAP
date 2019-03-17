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
#' matrix applying the transformations (see \code{whichtoinvert} in
#' \code{\link[ANTsRCore]{antsApplyTransforms}}).
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
#' @param templateRegMask antsImage or filename of the template mask that includes
#'               the skull but no face. Useful for improving the skull
#'               stripping process.
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
#'               saved at the specified path/prefix. The folder must exist or you
#'               will get an error. It is passed without modification to
#'               \code{antsRegistration}.
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
#'  \item\code{subImg} - subject\'s image in native space (after bias correction,
#'   denoising, skull stripping, etc.)
#'  \item\code{subLesion} - subject\'s lesion map in native space
#'  \item\code{subImgTemplate} - subject\'s image in template space
#'  \item\code{subLesionTemplate} - subject\'s lesion in template space
#'  \item\code{subRegMask} - registration mask in native space
#'  \item\code{templateImg} - the template used to register the subject
#'  \item\code{templateBrainMask} - the brain mask of the template image
#'  \item\code{subLesionTemplate} - the template mask with skull and no face
#'  \item\code{registration$inverse_subject2template} - transformation matrices subject to template
#'  \item\code{registration$forward_template2subject} - transformation matrices template to subject
#' }
#'
#' @examples
#' \dontrun{
#' anatomical = '/mnt/c/User/dp/Desktop/Subject1_anat.nii.gz'
#' lesion = '/mnt/c/User/dp/Desktop/Subject1_les.nii.gz'
#' newles = registerLesionToTemplate(anatomical, lesion,
#'         outprefix = '/mnt/c/User/dp/Desktop/Subj1onTemplate_')
#' }
#'
#' @author Dorian Pustina
#'

registerLesionToTemplate <- function(subImg, subLesion,
                                     templateImg = NA,
                                     templateBrainMask = NA,
                                     templateRegMask = NA,
                                     skullStrip=T,
                                     typeofTransform = 'SyNCC',
                                     outprefix = '',
                                     showInfo = T,
                                      ...) {


  toc = Sys.time()

  # load provided MNI template if user did not specify
  if (is.na(templateImg)) {
    if (showInfo) printInfo('Template undefined, using provided MNI152_2009c')
    mnipath = file.path(find.package('LESYMAP'), 'extdata', 'template','other_templates', 'MNI152_2009c')
    if (showInfo) printInfo(paste0('     ', mnipath), type='tail' )
    templateImg = antsImageRead(file.path(mnipath,'mni_icbm152_t1_tal_nlin_sym_09c.nii.gz'))
    templateBrainMask = antsImageRead(file.path(mnipath,'mni_icbm152_t1_tal_nlin_sym_09c_mask.nii.gz'))
    templateRegMask = antsImageRead(file.path(mnipath,'mni_icbm152_t1_tal_nlin_sym_09c_mask_skullnoface.nii.gz'))
    checkMask(templateImg, templateBrainMask)
  }


  # load images if filenames are passed
  if ( checkAntsInput(subImg) == 'antsFiles') {
    if (showInfo) printInfo('Loading subject\'s anatomical...')
    subImg = antsImageRead(subImg)
  }
  if ( checkAntsInput(subLesion) == 'antsFiles') {
    if (showInfo) printInfo('Loading subject\'s lesion file...')
    subLesion = antsImageRead(subLesion)
  }
  if ( checkAntsInput(templateImg) == 'antsFiles') {
    if (showInfo) printInfo('Loading template anatomical...')
    templateImg = antsImageRead(templateImg)
  }
  if ( !is.na(templateBrainMask) && checkAntsInput(templateBrainMask) == 'antsFiles') {
    if (showInfo) printInfo('Loading template brain mask...')
    templateBrainMask = antsImageRead(templateBrainMask)
  }
  if ( !is.na(templateRegMask) && checkAntsInput(templateRegMask) == 'antsFiles') {
    if (showInfo) printInfo('Loading template registration mask...')
    templateRegMask = antsImageRead(templateRegMask)
  }

  if ( is.na(templateBrainMask) & skullStrip) {
    if (showInfo) printInfo('templateBrainMask not specified. No skull stripping will be performed...')
    skullStrip = FALSE
  }



  # binarize mask inputs to be always 0-1
  subLesion = thresholdImage(subLesion, 0.1, Inf)
  if ( !is.na(templateBrainMask) ) templateBrainMask = thresholdImage(templateBrainMask, 0.1, Inf)
  if ( !is.na(templateRegMask) ) templateRegMask = thresholdImage(templateRegMask, 0.1, Inf)

  # make sure lesion and anatomical are in the same space
  if (showInfo) printInfo('Assuring lesion and antomical are in the same space...')
  checkMask(subImg, subLesion)


  # truncate outlier intensities, bias correct, and denoise
  if (showInfo) printInfo('Running bias correction on anatomical...')
  subImg = abpN4(subImg)
  if (showInfo) printInfo('Denoising the anatomical...')
  subImg = subImg %>% iMath('PeronaMalik', 10, 0.4)


  # here is the rationale
  # after skull stripping we take the mask and dilate it a bit
  # then we remove the lesion from it.
  # At the registration we use the brain-only image coupled with the
  # mask that is slightly larger. The reason for it is because we want
  # edge information (brain-air boundary) to improve the registration.
  # The approach is based on Nick Tustison's advice, see ANTsR issue:
  # https://github.com/ANTsX/ANTs/issues/483
  #
  if (skullStrip) {
    if (showInfo) printInfo('Skull-stripping subject\'s image...')
    temp = abpBrainExtraction(img = subImg, tem = templateImg,
                              temmask = templateBrainMask,
                              temregmask = templateRegMask,
                              regtype = 'SyNabp')

    subBrainMask = temp$bmask
    subBrainMask = thresholdImage(subBrainMask + subLesion, 0.5, Inf)
    subBrainMaskMD = subBrainMask %>% iMath('MD', 2)
    subRegMask = subBrainMaskMD - subLesion
    subImg = subImg * subBrainMask
    templateImg = templateImg * templateBrainMask

    if (nchar(outprefix) > 0)
      antsImageWrite(subBrainMask, filename = paste0(outprefix, '_subBrainMask.nii.gz'))
    if (nchar(outprefix) > 0)
      antsImageWrite(subImg, filename = paste0(outprefix, '_subBrainOnly.nii.gz'))

  } else {
    subRegMask = subImg * 0 + 1 - subLesion
  }


  # save registration mask
  if (nchar(outprefix) > 0)
    antsImageWrite(subRegMask, filename = paste0(outprefix, '_subRegMask.nii.gz'))


  # run registration
  # subject needs to be fixed and template moving
  # otherwise can't use mask with existing ANTsR
  if (showInfo) {
    printInfo(paste('Running template registration',
              ifelse(typeofTransform=='SyNCC','(expect >1 hour)',''),'...'))
  }
  reg = antsRegistration(fixed = subImg, moving = templateImg, mask = subRegMask,
                   outprefix = outprefix,
                   typeofTransform = typeofTransform)
  # save image if requested
  if (nchar(outprefix) > 0)
    antsImageWrite(reg$warpedfixout, filename = paste0(outprefix, '_subImgTemplate.nii.gz'))

  # bring lesion in template space
  if (showInfo) printInfo('Applying registration to lesion image...')

  subLesionTemplate = antsApplyTransforms(fixed = templateImg, moving = subLesion,
                                      transformlist = reg$invtransforms, whichtoinvert = c(1,0),
                                      interpolator = 'nearestNeighbor')
  # save lesion if requested
  if (nchar(outprefix) > 0)
    antsImageWrite(subLesionTemplate, filename = paste0(outprefix, '_subLesionTemplate.nii.gz'))

  if (showInfo) {
    printInfo(paste('Lesion size native:',
              round( sum(subLesion) * prod(antsGetSpacing(subLesion))/1000, 2),'ml'))
    printInfo(paste('Lesion size template:',
              round( sum(subLesion) * prod(antsGetSpacing(subLesionTemplate))/1000, 2),'ml'))
  }

  # prepare output list
  output = list()
  output$subImg = subImg
  output$subLesion = subLesion
  output$subImgTemplate = reg$warpedfixout
  output$subLesionTemplate = subLesionTemplate
  if (skullStrip) output$subBrainMask = subBrainMask
  output$subRegMask = subRegMask
  output$templateImg = templateImg
  output$templateBrainMask = templateBrainMask
  output$templateRegMask = templateRegMask
  output$registration = list()
  output$registration$inverse_subject2template = reg$invtransforms
  output$registration$forward_template2subject = reg$fwdtransforms

  tic = Sys.time()
  runtime = paste(round(as.double(difftime(tic,toc)),1), units(difftime(tic,toc)))
  if (showInfo) print(paste('Done!',runtime))

  return(output)

}
