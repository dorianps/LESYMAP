#' Support Vector Regression for symptom mapping
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' The SVR method is used. The function relies on the
#' \code{\link[e1071]{svm}} function of the e1071 package. The
#' analysis follows a similar logic found in the SVR-LSM code published
#' by \href{https://www.ncbi.nlm.nih.gov/pubmed/25044213}{Zhang (2015)}.
#' After a first run of SVM, p-values are established with a
#' permutation procedures as the number  of times weights are randomly
#' exceeded in permutations. The returned p-values are not corrected
#' for multiple comparisons.
#'
#' @param lesmat matrix of voxels (columns) and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param SVR.nperm (default=10,000) number of permutations to run to
#' estimate p-values. Note, these p-values are uncorrected for
#' multiple comparisons.
#' @param SVR.type (default='eps-regression') type of SVM to run, see
#' \code{\link[e1071]{svm}}.
#' @param SVR.kernel (default='radial') type of kernel to use, see
#' \code{\link[e1071]{svm}}
#' @param SVR.gamma (default=5) gamma value, see
#' \code{\link[e1071]{svm}}.
#' @param SVR.cost (default=30) cost value, see
#' \code{\link[e1071]{svm}}.
#' @param SVR.epsilon (default=0.1) epsilon value, see
#' \code{\link[e1071]{svm}}.
#' @param showInfo logical (default=TRUE) display messages
#' @param tstamp (default="\%H:\%M:\%S") timestamp format for messages
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{statistic} - vector of statistical values
#' \item\code{pvalue} - vector of pvalues
#' }
#'
#' @author Daniel Wiesen, Dorian Pustina
#'
#' @export
lsm_svr <- function(lesmat, behavior, mask,
                    SVR.nperm = 10000,
                    SVR.type = 'eps-regression',
                    SVR.kernel = 'radial',
                    SVR.gamma = 5,
                    SVR.cost = 30,
                    SVR.epsilon = 0.1,
                    showInfo = TRUE,
                    tstamp = "%H:%M:%S",
                    ...) {


  if (! 'e1071' %in% rownames(installed.packages())) stop('SVR-LSM requires e1071 package. Try installing with install.packages("e1071")')

  # library(e1071)

  # if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Write mask to file... ', "\n"))
  #
  # antsImageWrite(mask, filename = file.path(saveDir, 'mask.nii'))

  #get the idx to be able to create brain image after SVR-LSM
  # maskIdx = which(mask>0)
  # maskInfo = antsImageHeaderInfo(mask)
  # maskDim = maskInfo$dimensions
  #
  # #scaling behavioral
  # behavior = behavior*100/max(abs(behavior));

  # scale and center data
  behavior = scale(behavior, scale=T, center=T)
  lesmat = scale(lesmat, scale=T, center=T)

  # compute SVR
  if (showInfo) {
    cat(paste('\n       Calling SVR with:'))
    cat(paste('\n            SVM type =', SVR.type))
    cat(paste('\n            Kernel =', SVR.kernel))
    cat(paste('\n            Gamma =', SVR.gamma))
    cat(paste('\n            Cost =', SVR.cost))
    cat(paste('\n            Epsilon =', SVR.epsilon))
  }
  tic = Sys.time() # start time

  svr = e1071::svm(x = lesmat, y = behavior,
                   scale = FALSE, type = SVR.type,
                   kernel = SVR.kernel, gamma = SVR.gamma,
                   cost = SVR.cost, epsilon = SVR.epsilon,
                   na.action = na.fail)

  onerun = as.double(difftime(Sys.time(),tic, units = 'sec')) # total runtimme


  # get weights
  w = t(svr$coefs) %*% svr$SV

  # scale weights
  #' WHAT IS THE RATIONALE FOR SCALING BY 10?
  betaScale = 10/max(abs(w))
  statistic = as.vector(w*betaScale)



  # run permutations to establish p-values for weights (uncorrected)
  if (showInfo) {
    toc = tic + (SVR.nperm * onerun)
    expect = paste(round(as.double(difftime(toc,tic)),1), units(difftime(toc,tic)))
    cat(paste('\n       SVR', SVR.nperm, 'permutations, expected run =', expect))
  }

  exceed = rep(1, length(statistic))
  posindx = statistic >=0
  negindx = statistic < 0
  for (i in 1:SVR.nperm) {

    svr = e1071::svm(x = lesmat, y = sample(behavior, replace=FALSE),
                     scale = FALSE, type = SVR.type,
                     kernel = SVR.kernel, gamma = SVR.gamma,
                     cost = SVR.cost, epsilon = SVR.epsilon,
                     na.action = na.fail)


    # get weights
    w = t(svr$coefs) %*% svr$SV

    # scale weights
    #' WHAT IS THE RATIONALE FOR SCALING BY 10?
    betaScale = 10/max(abs(w))
    thisstat = as.vector(w*betaScale);

    exceed[posindx] = exceed[posindx] + (thisstat[posindx] >= statistic[posindx])
    exceed[negindx] = exceed[negindx] + (thisstat[negindx] <= statistic[negindx])

  }

  pvalue = exceed / (SVR.nperm + 1)

  output = list()
  output$statistic = statistic
  output$pvalue = pvalue

  return(output)







  #' DP: I AM LOST BELOW THIS POINT
  #' THE FUNCTION SHOULD SIMPLY RETURN THE BETA WEIGHTS
  #' AND EVENTUAL P VALUES
  #' LESYMAP TAKES CARE OF PUTTING VALUES BACK INTO AN IMAGE


  # #preallocate memory and get indices
  # betaMap = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # posBetaMap = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # negBetaMap = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # betaMapAbs = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  #
  # #create betaMaps
  # betaMap[maskIdx] = beta_weights
  # posIdx = which(betaMap>=0)
  # negIdx = which(betaMap<0)
  # pMapPosIdx = which(betaMap>0)
  # posBetaMap[posIdx] = betaMap[posIdx]
  # negBetaMap[negIdx] = betaMap[negIdx]
  # betaMapAbs = abs(betaMap)
  #
  # #create beta images to write out
  # betaImage = makeImage(c(maskDim[1],maskDim[2],maskDim[3]), voxval= betaMapAbs)
  # posBetaImage = makeImage(c(maskDim[1],maskDim[2],maskDim[3]), voxval= posBetaMap)
  # negBetaImage = makeImage(c(maskDim[1],maskDim[2],maskDim[3]), voxval= negBetaMap)
  # antsCopyImageInfo(mask, betaImage)
  # antsCopyImageInfo(mask, posBetaImage)
  # antsCopyImageInfo(mask, negBetaImage)
  #
  # if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Write Out Beta Map...', "\n"))
  #
  # #write out beta images
  # antsImageWrite(betaImage, filename = file.path(saveDir, 'beta_map_abs.nii'))
  # antsImageWrite(posBetaImage, filename = file.path(saveDir, 'beta_map_pos.nii'))
  # antsImageWrite(negBetaImage, filename = file.path(saveDir, 'beta_map_neg.nii'))

  #### run permutations

  # if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Run Permutations...', "\n"))
  #
  # #preallocate memory
  # map_count = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # map_count_pos = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # map_count_neg = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # ori_beta_val = betaMap
  #
  # #start permutation
  #
  # for(PermIdx in 1:nperm) {
  #
  #   maps = run_beta_PMU(lesmat, behavior, ori_beta_val, maskIdx, posIdx, negIdx, maskDim, betaScale)
  #
  #   map_count_tmp = maps$map_count_tmp
  #   map_count_pos_tmp = maps$map_count_pos_tmp
  #   map_count_neg_tmp = maps$map_count_neg_tmp
  #
  #   map_count = map_count + map_count_tmp   #update the counter with current permutation for producing p_map
  #   map_count_pos = map_count_pos + map_count_pos_tmp #update the counter with current permutation for producing p_map
  #   map_count_neg = map_count_neg + map_count_neg_tmp #update the counter with current permutation for producing p_map
  #
  # }
  #
  # if (showInfo) cat(paste(format(Sys.time(), tstamp) , 'Generate and Write Out untresholded p_maps...', "\n"))
  #
  # #generate p_maps
  # p_map = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # p_map = map_count/(nperm+1);
  # p_map_pos = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # p_map_pos = map_count_pos/(nperm+1);
  # p_map_neg = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # p_map_neg = map_count_neg/(nperm+1);
  #
  # #flip to 1-p to be able to visualize in mricron, create image and write out to file
  # tmp_out_img = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # tmp_out_img[maskIdx] = 1-p_map[maskIdx];
  # tmp_toWrite = makeImage(c(maskDim[1],maskDim[2],maskDim[3]), voxval= tmp_out_img)
  # antsCopyImageInfo(mask, tmp_toWrite)
  # antsImageWrite(tmp_toWrite, filename = file.path(saveDir, 'untresholded_p_map_inv.nii'))
  #
  # tmp_out_img = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # tmp_out_img[pMapPosIdx] = 1-p_map_pos[pMapPosIdx];
  # tmp_toWrite = makeImage(c(maskDim[1],maskDim[2],maskDim[3]), voxval= tmp_out_img)
  # antsCopyImageInfo(mask, tmp_toWrite)
  # antsImageWrite(tmp_toWrite, filename = file.path(saveDir, 'untresholded_p_map_pos_inv.nii'))
  #
  # tmp_out_img = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
  # tmp_out_img[negIdx] = 1-p_map_neg[negIdx];
  # tmp_toWrite = makeImage(c(maskDim[1],maskDim[2],maskDim[3]), voxval= tmp_out_img)
  # antsCopyImageInfo(mask, tmp_toWrite)
  # antsImageWrite(tmp_toWrite, filename = file.path(saveDir, 'untresholded_p_map_neg_inv.nii'))


}

