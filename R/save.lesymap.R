#' Save the output of lesymap.
#'
#' Function to save the output of lesymap.
#'
#' @param lsm object obtained with lesymap()
#' @param saveDir folder to save to, will be
#' created if it doesn't exist.
#' @param infoFile (default='Info.txt') what should
#' be the filename of the file with information.
#' @param template (default=NA) an antsImage to overlay
#' the results to. If the template is provided,
#' results will be plotted and saved as image.
#' @param saveTemplate (default=FALSE) should the
#' template image also be saved? Useful when
#' passing the results to a colleague.
#' @param plot.alpha see plot.antsImage
#' @param plot.axis see plot.antsImage
#' @param plot.quality see plot.antsImage
#' @param ... other arguments to use for plot().
#'
#' @return Nothing is returned.
#' Files saved include both images and
#' a descriptive file with a lot of
#' information about the lesymap run.
#'
#' @author Dorian Pustina
#'
#' @export
save.lesymap <- function(lsm, saveDir, infoFile='Info.txt', template=NA, saveTemplate=F,
                          plot.alpha=0.8, plot.axis=3, plot.quality=8, ...) {

  if ('callinfo' %in% names(lsm)) callinfo = lsm$callinfo
  else callinfo=list()

  # prepare output folder
  if (! file.exists(saveDir) ) {
    boolcreate = dir.create(saveDir, showWarnings = F, recursive = T)
    if (!boolcreate) {
      warning(paste('Cannot create saving folder:', saveDir, '\nAborting save.'))
      return()
    }
  }

  # do we have a template?
  hastemplate = F
  if (!is.na(template)) {
    if (checkAntsInput(template) == 'antsFiles') template = antsImageRead(template)
    if (!checkImageList(list(template, lsm$stat.img), showError = F, binaryCheck = F)) {
      warning('Ignoring template: not in same space with lesymap output.')
    } else {
      hastemplate = T
    }
  }


  # save images
  imgindx =  grep('.img$', names(lsm) )
  for (indx in imgindx) {

    if (class(lsm[[indx]]) != 'antsImage') next # smth wrong if this is not image

    thisname = names(lsm)[indx]

    antsImageWrite(lsm[[indx]], filename = file.path(saveDir, paste0(thisname, '-.nii.gz') ))

    if (thisname == 'stat.img' & hastemplate) {
      outname = file.path(saveDir, paste0(thisname, '.png') )
      plot.antsImage(template, lsm[[indx]], alpha=plot.alpha, axis=plot.axis,
                     window.overlay=range(lsm[[indx]]), useAbsoluteScale=F,
                     outname=outname, quality=plot.quality, ...)
    }
  }
  if (hastemplate & saveTemplate) antsImageWrite(template, filename = file.path(saveDir, 'template.nii.gz'))


  # save info
  outfile =  file.path(saveDir, infoFile)
  file.create(outfile)

  # write lesymap call parameters
  if (length(callinfo)>0) {
    for (i in 1:length(callinfo)) {
      line = paste0(names(callinfo)[i], ': ', callinfo[[i]])
      write(line, outfile, append = T)
    }
  }

  # write patch information
  if ('patchinfo' %in% names(lsm)) {
    line = paste0('Total voxels in mask: ', lsm$patchinfo$nvoxels)
    write(line, outfile, append = T)
    line = paste0('Unique voxel patches: ', lsm$patchinfo$npatches)
    write(line, outfile, append = T)
    line = paste0('Average patch size: ', round(mean(lsm$patchinfo$patchvoxels), 1), ' voxels')
    write(line, outfile, append = T)
    line = paste0('Unique single voxels: ', sum(lsm$patchinfo$patchvoxels==1) )
    write(line, outfile, append = T)
  }

  line = paste('Range statistic:', paste( round(range(lsm$stat.img),2) , collapse=' '))
  write(line, outfile, append = T)
  line = paste('Suprathreshold voxels:', sum(lsm$stat.img!=0))
  write(line, outfile, append = T)

  clust.stat = labelStats(lsm$stat.img, labelClusters(abs(lsm$stat.img), minThresh = .Machine$double.eps, maxThresh = Inf))
  clust.stat = clust.stat[clust.stat$LabelValue!=0,]
  line = paste('Number of Clusters:', nrow(clust.stat))
  write(line, outfile, append = T)
  line = paste('Cluster sizes (voxels):', paste(clust.stat$Count, collapse=' '))
  write(line, outfile, append = T)
  line = paste('Cluster sizes (mm3):', paste(clust.stat$Volume, collapse=' '))
  write(line, outfile, append = T)

  # write other occasional outputs
  # must be a single scalar value in lsm to be written
  for (i in 1:length(lsm)) {
    if (is.numeric(lsm[[i]]) & length(lsm[[i]]) <= 2) {
      line = paste0(names(lsm)[i], ': ', paste(lsm[[i]], collapse=' '))
      write(line, outfile, append = T)
    }
    if (is.character(lsm[[i]]) & length(lsm[[i]]) <= 2) {
      line = paste0(names(lsm)[i], ': ', paste(lsm[[i]], collapse=' '))
      write(line, outfile, append = T)
    }
  }

}
