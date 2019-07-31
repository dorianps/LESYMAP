#' @title Save the output of lesymap.
#'
#' @description
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
#' @param savePatchImages (default=TRUE) should the patch images
#' be saved
#' @param plot.alpha see plot.antsImage
#' @param plot.axis see plot.antsImage
#' @param plot.quality see plot.antsImage
#' @param outputLogFile (default='outputLog.txt')
#' the filename to save the console output
#' @param ... other arguments to use for plot().
#'
#' @return Nothing is returned.
#' Files saved include resulting maps and
#' a descriptive file with a lot of
#' information about the lesymap run.
#'
#' @author Dorian Pustina
#'
#' @export
#' @importFrom graphics plot
save.lesymap <- function(lsm, saveDir, infoFile='Info.txt', template=NA, saveTemplate=F,
                          savePatchImages=T, plot.alpha=0.8, plot.axis=3, plot.quality=8,
                         outputLogFile = 'outputLog.txt', ...) {

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
  if (!is.na(c(template))) {
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
    savename = gsub('.img$', '_img', thisname)

    antsImageWrite(lsm[[indx]], filename = file.path(saveDir, paste0(savename, '.nii.gz') ))

    if (thisname == 'stat.img' & hastemplate) {
      outname = file.path(saveDir, paste0(thisname, '.png') )
      # plot.antsImage(template, lsm[[indx]], alpha=plot.alpha, axis=plot.axis,
      #                window.overlay=range(lsm[[indx]]), useAbsoluteScale=F,
      #                outname=outname, quality=plot.quality, ...)
      plot(template, lsm[[indx]], alpha=plot.alpha, axis=plot.axis,
                     window.overlay=range(lsm[[indx]]), useAbsoluteScale=FALSE,
                     outname = outname, quality=plot.quality, ...)
    }
  }
  if (hastemplate & saveTemplate) antsImageWrite(template, filename = file.path(saveDir, 'template.nii.gz'))


  # patchinfo image save
  if (savePatchImages & 'patchinfo' %in% names(lsm)) {
    imgindx =  grep('^patchimg', names(lsm$patchinfo) )
    for (indx in imgindx) {
      if (class(lsm[[indx]]) != 'antsImage') next # smth wrong if this is not image
      thisname = names(lsm$patchinfo)[indx]
      antsImageWrite(lsm$patchinfo[[indx]], filename = file.path(saveDir, paste0(thisname, '.nii.gz') ))
    }
  }



  # save info
  outfile =  file.path(saveDir, infoFile)
  file.create(outfile)

  # write lesymap call parameters
  if (length(callinfo)>0) {
    for (i in 1:length(callinfo)) {
      if (is.numeric(callinfo[[i]])) {
        thisval = paste(callinfo[[i]], collapse=' ')
      } else if (typeof(callinfo[[i]]) == 'language') { # special case for negative numbers
        thisval = paste(callinfo[[i]], collapse='')
      } else {
        thisval = callinfo[[i]]
      }
      line = paste0(names(callinfo)[i], ': ', thisval)
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

  if ("rawWeights.img" %in% names(lsm)) {
    line = paste('Range rawWeights:', paste( range(lsm$rawWeights.img) , collapse=' '))
    write(line, outfile, append = T)
  }

  if (sum(lsm$stat.img!=0) > 0) {
    clust.stat = labelStats(lsm$stat.img, labelClusters(abs(lsm$stat.img), minThresh = .Machine$double.eps, maxThresh = Inf))
    clust.stat = clust.stat[clust.stat$LabelValue!=0,]
    line = paste('Number of Clusters:', nrow(clust.stat))
    write(line, outfile, append = T)
    line = paste('Cluster sizes (voxels):', paste(clust.stat$Count, collapse=' '))
    write(line, outfile, append = T)
    line = paste('Cluster sizes (mm3):', paste(clust.stat$Volume, collapse=' '))
    write(line, outfile, append = T)
  }

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
    if (names(lsm)[i] == "sccan.ccasummary") {
      line = paste0(names(lsm)[i], ': ', paste(as.numeric(lsm[[i]])[1], collapse=' '))
      write(line, outfile, append = T)
    }
  }

  # write printOutput
  if ('outputLog' %in% names(lsm)) {
    outlog =  file.path(saveDir, outputLogFile)
    file.create(outlog)
    for (i in 1:length(lsm$outputLog)) write(lsm$outputLog[i], outlog, append = T)
  }

}
