#' Lesion to Symptom Mapping
#'
#' Lesymap uses univariate and multivariate methods to map
#' functional regions of the brain that, when lesioned,
#' cause specific cognitive deficits. It requires is
#' a set of binary lesion maps (nifti files) in template
#' space and the vector of corresponding behavioral scores.
#' Note, lesions must be already registered in template space,
#' use the built-in function
#' \code{\link[=registerLesionToTemplate]{registerLesionToTemplate}}
#' or other ANTs tools to register lesions. Lesymap will check that
#' lesion maps are in the same space before running. For traditional
#' mass-univariate analyses (i.e., BMfast, ttest, etc.),
#' voxels with identical lesion patterns are grouped together in
#' unique patches. Patch-based analysis decreases the number of
#' multiple comparisons and speeds up the analyses. Multivariate
#' analysis are performed using an optimized version of sparse
#' canonical correlations (SCCAN, \code{method='sccan'})or support
#' vector regression (\code{method='svr'}).
#'
#' @param lesions.list list of antsImages, or a vector of
#' filenames, or a single antsImage with 4 dimensions.
#'
#' @param behavior vector of behavioral scores or filename
#' pointing to a file with a single column of numbers.
#'
#' @param method  what analysis method to use to run, one of
#' 'BM', 'BMfast', 'ttest', 'welch', 'regres', 'regresfast',
#' 'regresPerm', 'sccan' (default) or 'svr',.
#'
#'        \code{\link[=lsm_BM]{BM}} - Brunner-Munzel non parametric test, also
#'          called the Generalized Wilcoxon Test. The BM test is the
#'          same test used in the npm/Mricron software
#'          (see \href{https://www.ncbi.nlm.nih.gov/pubmed/17583985}{Rorden (2007)}).
#'          This method is slow, use 'BMfast" for a compiled faster analysis.
#'
#'        \code{\link[=lsm_BMfast]{BMfast}} - ultrafast Brunner-Munzel with compiled
#'          code. BMfast can be
#'          combined with \code{multipleComparison='FWERperm'}
#'          to perform permutation based thresholding in a short time.
#'
#'        \code{\link[=lsm_ttest]{ttest}} - Regular single tailed t-test.
#'          Variances of groups are assumed to be equal at each voxel, which may not
#'          be true. This is the test used in the voxbo
#'          software. Relies on \code{t.test} function in R. It is assumed
#'          that 0 voxels are healthy, i.e., higher behavioral scores.
#'          See the \code{alternative} parameter for inverted cases.
#'          (see \href{https://www.ncbi.nlm.nih.gov/pubmed/12704393}{Bates (2003)}).
#'
#'        \code{\link[=lsm_ttest]{welch}} - t-test that does not assume equal
#'          variance between groups. Relies on \code{t.test} function in R.
#'
#'        \code{\link[=lsm_regres]{regres}} - linear model between voxel values
#'          and behavior. Uses the \code{lm} function in R. This is equivalent to a
#'          t-test, but is useful when voxel values are continuous. This method is
#'          R-based and slow, for faster analysis and to add covariates
#'          use the \code{"regresfast"} method compiled in LESYMAP.
#'
#'        \code{\link[=lsm_regresfast]{regresfast}} - fast linear regressions with
#'          compiled code. This method allows setting covariates. If covariates are
#'          specified the effect of each voxel will be estimated with the formula:\cr
#'          \code{behavior ~ voxel + covar1 + covar2 + ...}\cr
#'          The effect of covariates at each voxel is established with
#'          the Freedman-Lane method
#'          (see \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4010955/}{Winkler (2014)}).
#'          This LESYMAP method allows multiple comparison correction with permutation
#'          based methods \code{"FWERperm"} and \code{"clusterPerm"}.
#'
#'        \code{\link[=lsm_regresPerm]{regresPerm}} - linear model between voxel values
#'          and behavior. The p-value of each individual voxel is established by
#'          permuting voxel values. The lmPerm package is used for this purpose.
#'          Note, these permutations do not correct for
#'          multiple comparisons, they only establish voxel-wise p-values, a correction
#'          for multiple comparisons is still required.
#'
#'        \code{\link[=lsm_chisq]{chisq}} - chi-square test between voxel values and behavior.
#'          The method is used when behavioral scores are binary (i.e. presence of absence
#'          of deficit). Relies on the \code{\link[stats]{chisq.test}}
#'          R function. By default this method corrects individual voxel p-values with the
#'          Yates method (the same approach offered in the Voxbo software).
#'
#'        \code{\link[=lsm_chisq]{chisqPerm}} - chi-square tests. P-values are established
#'          through permutation tests instead of regular statistics.
#'          Relies on the \code{\link[stats]{chisq.test}} R function.
#'
#'        \code{\link[=lsm_sccan]{sccan}} - sparse canonical correlations. Multivariate
#'          method that considers all voxels at once by searching for voxel weights
#'          that collectively explain behavioral variance. For our purposes, this
#'          method can be considered a sparse regression technice. By default,
#'          lesymap will run a lengthy procedure to determine the optimal
#'          sparseness value (how extensive the results should be). You
#'          can set \code{optimizeSparseness=FALSE} if you want to skip this
#'          optimization. The search for optimal sparsness provides
#'          a cross-validated correlation measure that shows how well the
#'          sparseness value can predict new patients. If this predictive correlation
#'          is below significance (i.e., below \code{pThreshold}), the entire solution
#'          will be dropped and LESYMAP will return an empty statistical map.
#'          If the predictive correlation is significant, LESYMAP will return
#'          a statistical image with normalized weights between -1 and 1. The raw SCCAN
#'          weights image is returned in \code{rawWeights.img} as well. LESYMAP scales
#'          and centers both lesion and behavior data before running SCCAN
#'          (hardcoded in \link{lsm_sccan}). See more details in
#'          \href{https://www.ncbi.nlm.nih.gov/pubmed/28882479}{Pustina (2018)}
#'
#'        \code{\link[=lsm_svr]{svr}} - support vector regression. Multivariate
#'          method that considers all voxels at once by searching for voxel weights.
#'          To establish p-values, a number of permutations are needed as set with
#'          \code{SVR.nperm}. The current implementation of SVR in LESYMAP is not
#'          parallelized and takes many hours to finish all permutations. The SVR
#'          method is initially described in
#'          \href{https://www.ncbi.nlm.nih.gov/pubmed/25044213}{Zhang (2014)},
#'          LESYMAP uses a contribution by the
#'          \href{https://www.ncbi.nlm.nih.gov/pubmed/30549154}{Tuebingen group}.
#'
#' @param minSubjectPerVoxel (default='10\%') remove voxels/patches with lesions
#'  in less than X subjects. Value can
#'  be speficifed as percentage ('10\%')
#'  or exact number of subjects (10).
#'
#' @param correctByLesSize whether to correct for lesion size in the analysis.
#'  Options are "none", "voxel", "behavior", "both":
#'  \itemize{
#'  \item\code{"none"}: (default) no correction
#'  \item\code{"voxel"}: divide voxel values by 1/sqrt(lesionsize).
#'     This is the method used in
#'     \href{https://www.ncbi.nlm.nih.gov/pubmed/25879574}{Mirman (2015)} and
#'     \href{https://www.ncbi.nlm.nih.gov/pubmed/25044213}{Zhang (2014)}.
#'     This correction works only with
#'     'regres' methods. Two sample comparisons
#'     (t-tests and Brunner-Munzel) use binary voxels
#'     and will ignore this correction.
#'  \item\code{"behavior"}: residualize behavioral scores by removing
#'     the effect of lesion size. This works on all methods,
#'     but is more agressive on results.
#'  \item\code{"both"}: both voxel and behavior residualized.
#'     }
#'
#' @param multipleComparison (default='fdr') method to adjust p-values.
#'  Standard methods include \code{"holm"}, \code{"hochberg"},
#'  \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}.
#'  (see \code{\link[stats]{p.adjust}} ) \cr
#'  Permutation methods include: \cr
#'  \code{"FWERperm"} (permutation based family-wise threshold) is enabled
#'    with methods 'BMfast' and 'regresfast'. In this
#'    case, many analysis are run with permuted behavioral
#'    scores, and the peak score is recorded each time (see
#'    Winkler 2014). The optimal threshold is established at
#'    95th percentile of this distribution (or whatever pThreshold
#'    you choose). You can choose to use as reference another
#'    voxel lower in the ranks by specifying another `v` value
#'    (i.e., lesymap(..., v=10) will record the 10th highest voxel). \cr
#'  \code{"clusterPerm"} (permutation based cluster correction) is enabled
#'    for 'regresfast'. It records the maximal
#'    cluster size from many random permutations of the behavior
#'    score and sets a cluster threshold based on that distribution.
#'    You must select pThreshold (voxel-wise, default=0.05) and
#'    clusterPermThreshold (cluster-wise, default 0.05) to achieve
#'    optimal thresholding with this approach.
#'
#' @param pThreshold (default=0.05) threshold statistics at this p-value
#'  (after corrections or permutations)
#'
#' @param mask (default=NA) binary image to select the area
#'  where analysis will be performed. If
#'  not provided will be computed automatically
#'  by thresholding the average lesion map at
#'  \code{minSubjectPerVoxel}.
#'
#' @param patchinfo (default=NA) an object obtained with the
#'  \code{getUniqueLesionPatches} function. Useful to set if your
#'  run repetitive analysese and want toavoid the computation of
#'  patches each time. You can also pass the patchinfo object
#'  obtained from a previous analysis.
#'
#' @param nperm (default=1000) number of permutations to run when
#'   necessary. This is used mostly for univariate analyses, while
#'   multivariate methods have their own permutation arguments.
#'   Check the documentation of each method to know more.
#'
#' @param noPatch logical (default=FALSE), if True avoids using patch information
#'  and will analyze all voxels. It will take longer and results will be
#'  worse due to more multiple comparison corrections. This argument
#'  is ignored when performing multivariate analyses, SCCAN or SVR, for
#'  which all voxels are always used.
#'
#' @param flipSign logical (default=FALSE), invert the sign in the statistics image.
#'
#' @param binaryCheck logical (default=FALSE), make sure the lesion matrix is 0/1.
#'   This will help if lesion maps are drawn in MRIcron or other software which
#'   label lesioned voxel with value 255.
#'
#' @param showInfo logical (default=TRUE), display time-stamped info messages
#'
#' @param saveDir (default=NA) save results in the specified folder.
#'
#' @param ... arguments that will be passed down to other functions
#' (i.e., sparsness=0.045)
#'
#'
#' @details
#' Several other parameters can be specified to lesymap()
#' which will be passed to other called fuctions. Here are
#' some examples:
#'
#' \code{permuteNthreshold}  - (default=9) for Brunner-Munzel tests
#'    in \code{method='BMfast'} or \code{method='BM'}.
#'    Voxels lesioned in less than this number
#'    of subjects will undergo permutation-based
#'    p-value estimation. Useful because the BM test
#'    is not valid when comparing groups with N < 9.
#'    Issue described in:
#'    \href{https://www.ncbi.nlm.nih.gov/pubmed/19766664}{Medina (2010)}
#'
#' \code{clusterPermThreshold} - threshold used to find the optimal cluster size when
#'    using \code{multipleComparison='clusterPerm'}.
#'
#' \code{alternative} - (default='greater') for two sample tests (ttests and BM).
#'    By default LESYMAP computes single tailed p-values assuming
#'    that non-lesioned 0 voxels have higher behavioral scores.
#'    You can specify the opposite relationship with alternative='less'
#'    or compute two tailed p-values with alternative='two.sided'.
#'
#' \code{covariates} - (default=NA) enabled for method = 'regresfast'.
#'    This will allow to model the effect of each voxel
#'    in the context of other covariates, i.e., formula
#'    "behavior ~ voxel + covar1 + covar2 + ...". \cr
#'    I,.e., lesymap(lesions,behavior, method='regresfast',
#'            covariates=cbind(lesionsize, age)).\cr
#'    If you choose permutation based thresholding with covariates, lesymap
#'    will use the Freedman-Lane method for extracting the unique effect of
#'    each voxel (see Winkler 2014, Freedman 1983)
#'
#' \code{template} - antsImage or filename used for plotting the results if a saving
#'    directory is specified (see \code{saveDir})
#'
#' \code{v} - (default=1) which voxel to record for permutation based thresholding.
#'    Normally the peak voxel is used (1), but other voxels can be recorded.
#'    See \href{https://www.ncbi.nlm.nih.gov/pubmed/28847712}{Mirman (2017)}
#'    for this approach.
#'
#'
#'
#' @return
#' The following objects are typically found in the returned list:
#' \itemize{
#' \item\code{stat.img}   - statistical map
#' \item\code{rawWeights.img}   - (optional) raw SCCAN weights
#' \item\code{pval.img}   - (optional) p-values map
#' \item\code{zmap.img}   - (optional) zscore map
#' \item\code{mask.img}    - mask used for the analyses
#' \item\code{average.img} - map of all lesions averaged,
#'            produced only if no mask is defined.
#' \item\code{callinfo} - list of details of how you called lesymap
#' \item\code{outputLog} - terminal output in a character variable
#' \item\code{perm.vector} - (optional) the values obtained from each permutation
#' \item\code{perm.clusterThreshold} - (optional) threshold computed
#'            for cluster thresholding
#' \item\code{perm.FWERthresh} - (optional) threshold computed for FWERperm
#'            thresholding
#' \item\code{patchinfo}   - list of variables describing patch information:
#'     \itemize{
#'     \item\code{patchimg} - antsImage with the patch number each voxels belongs to
#'     \item\code{patchimg.samples} - antsImage mask with a single voxel per patch
#'     \item\code{patchimg.size} - antsImage with the patch size at each voxel
#'     \item\code{patchimg.mask} -  the mask within which the function will look for patches
#'     \item\code{npatches} - number of unique patches in the image
#'     \item\code{nvoxels} - total number of lesioned voxels in mask
#'     \item\code{patchvoxels} - vector of voxel count for each patch
#'     \item\code{patchvolumes} - vector of volume size for each patch
#'     \item\code{patchmatrix} - the lesional matrix, ready for use in analyses.
#'            Matrix has size NxP (N=number of subjects, P=number of
#'            patches)
#'            }
#' }
#'
#' @examples
#' lesydata = file.path(find.package('LESYMAP'),'extdata')
#' filenames = Sys.glob(file.path(lesydata, 'lesions', 'Subject*.nii.gz'))
#' behavior = Sys.glob(file.path(lesydata, 'behavior', 'behavior.txt'))
#' template = antsImageRead(
#'  Sys.glob(file.path(lesydata, 'template', 'ch2.nii.gz')))
#' lsm = lesymap(filenames, behavior, method = 'BMfast')
#' plot(template, lsm$stat.img, window.overlay = range(lsm$stat.img))
#'
#' \dontrun{
#' # Same analysis with SCCAN
#' lsm = lesymap(filenames, behavior, method = 'sccan',
#' sparseness=0.045, validateSparseness=FALSE)
#' plot(template, lsm$stat.img, window.overlay = range(lsm$stat.img))
#' save.lesymap(lsm, saveDir='/home/dp/Desktop/SCCANresults')
#' }
#'
#' @author Dorian Pustina
#'
#' @export
#' @useDynLib LESYMAP
#' @importFrom Rcpp sourceCpp
#' @import ANTsR
#' @import ANTsRCore
#' @importFrom stats chisq.test cor dt lm optimize na.fail
#' @importFrom stats p.adjust p.adjust.methods pt qchisq qnorm
#' @importFrom stats quantile residuals runif t.test
#' @importFrom utils find getSrcDirectory installed.packages
#' @importFrom utils packageVersion read.table capture.output



lesymap <- function(lesions.list, behavior,
                    mask=NA,
                    patchinfo=NA,
                    method='sccan',
                    correctByLesSize='none',
                    multipleComparison='fdr',
                    pThreshold=0.05,
                    flipSign=F,
                    minSubjectPerVoxel = '10%',
                    nperm=1000,
                    saveDir=NA,
                    binaryCheck=FALSE,
                    noPatch=FALSE, showInfo=TRUE,
                    ...) {
  ver = as.character(packageVersion('LESYMAP'))
  toc = Sys.time()

  # start capturing window output to save later
  outputLog = capture.output({

    if (showInfo) printInfo(paste('Running LESYMAP', ver))
    if (showInfo) printInfo('Checking a few things...')

    # the full array of methods and other arguments accepted by lesymap
    acceptedArgs = list()
    acceptedArgs$BM = list(multipleComparison = c(p.adjust.methods),
                           correctByLesSize = c('none', 'behavior') )
    acceptedArgs$BMfast = list(multipleComparison = c(p.adjust.methods, 'FWERperm'),
                               correctByLesSize = c('none', 'behavior') )
    acceptedArgs$ttest = list(multipleComparison = c(p.adjust.methods),
                              correctByLesSize = c('none', 'behavior') )
    acceptedArgs$welch = list(multipleComparison = c(p.adjust.methods),
                              correctByLesSize = c('none', 'behavior') )
    acceptedArgs$regres = list(multipleComparison = c(p.adjust.methods),
                               correctByLesSize = c('none', 'behavior', 'voxel', 'both') )
    acceptedArgs$regresfast = list(multipleComparison = c(p.adjust.methods, 'FWERperm', 'clusterPerm'),
                                   correctByLesSize = c('none', 'behavior', 'voxel', 'both') )
    acceptedArgs$regresperm = list(multipleComparison = c(p.adjust.methods),
                                   correctByLesSize = c('none', 'behavior', 'voxel', 'both') )
    acceptedArgs$chisq = list(multipleComparison = c(p.adjust.methods),
                              correctByLesSize = c('none') )
    acceptedArgs$chisqperm = list(multipleComparison = c(p.adjust.methods),
                                  correctByLesSize = c('none') )
    acceptedArgs$sccan = list(multipleComparison = c(p.adjust.methods),
                              correctByLesSize = c('none', 'behavior', 'voxel', 'both') )

    # now check the inputs
    methodID = match( tolower(method), tolower(names(acceptedArgs)) )
    if (is.na(methodID)) {
      stop(paste0('Unrecognized method: ', method, '. Available methods: ',
                  paste(names(acceptedArgs), collapse=', ')
                  ))
    }
    method = names(acceptedArgs)[methodID] # ignore case, use our internal
    for (b in 1:length(acceptedArgs[[methodID]])) {
      thisAcceptedArg = names(acceptedArgs[[methodID]][b])
      thisAcceptedValues = acceptedArgs[[methodID]][[b]]
      thisArgValue = eval(parse(text=thisAcceptedArg))
      argID = match(tolower(thisArgValue), tolower(thisAcceptedValues))

      if (is.na(argID)) {
        stop(paste0('Method ', method, ' - Unrecognized value for ',
                    thisAcceptedArg, '. Available options: ',
                    paste(thisAcceptedValues, collapse=', ')
                    ))
      }

      # ignore case for user, use our internal option
      eval(parse(text=paste0(thisAcceptedArg, '= thisAcceptedValues[argID]')))
    }
    # finished checking inputs

    # load behavior if filename
    if (checkAntsInput(behavior) == 'antsFiles') {
      if (showInfo) printInfo('Loading behavioral data...', type='head')
      behavior = read.table(behavior, header = FALSE)$V1
      if (showInfo) printInfo(paste(length(behavior), 'scores found.'), type='tail')
    }

    inputtype = checkAntsInput(lesions.list, checkHeaders = TRUE)

    # make sure input is a list or filenames
    if ( !(inputtype %in% c('antsImageList', 'antsFiles', 'antsImage')) ) {
      stop('Unrecognized input: lesions.list should be a 4D antsImage, a list of antsImages, or a vector of filenames')
    }

    ##############
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

    #########
    # check lesions and behavior have same length
    if (length(lesions.list) != length(behavior)) stop('Different lengths between lesions and behavior vector.')

    ########
    # check behavior is binary if needed
    if (method %in% c('chisq', 'chisqPerm')) {
      if (length(unique(behavior)) != 2) stop(paste0('The method "', method, '" requries binary behavioral scores.'))
    }

    ########
    # special case for SCCAN
    if (method %in% c('sccan', 'sccanRaw') ) {
      if (showInfo) printInfo('SCCAN method: ignoring patch, nperm, and multiple comparison...')
      multipleComparison = 'none'
      noPatch = TRUE
      nperm=0
    }

    ########
    # special case for SVR
    if (method %in% c('svr') ) {
      if (showInfo) printInfo('SVR method: ignoring nperm, use SVR.nperm instead...')
      nperm=0
    }


    ## CHECK AND PREPARE MASK
    writeavgles = F
    if (!is.na(mask)) { # USER DEFINED

      if (showInfo) printInfo('Using predefined mask...')

      # check mask and lesions have same headers
      checkMask(lesions.list, mask) # first, if mask is a file
      if (checkAntsInput(mask) == 'antsFiles') { # load the mask if file
        mask=antsImageRead(mask)
        if (class(mask) != 'antsImage') stop('Mask must be of class antsImage')
        checkMask(lesions.list, mask) # check again after loading in memory
      }

      # check mask is binary
      if (! checkImageList(list(mask), binaryCheck = T, showError = F))
        stop('Mask is not binary. Must have only 0 and 1 values')

    } else if (!is.na(patchinfo[1])) { # GET IT FROM patchinfo

      if (showInfo) printInfo('Using mask passed with patchinfo...')
      mask = patchinfo$patchimg.mask
      checkMask(lesions.list, mask)

    } else { # DEFINE THE MASK
      if (showInfo) printInfo(paste('Searching voxels lesioned in >=',minSubjectPerVoxel,'subjects...'), type='head')
      # compute thresholdPercent based on minSubjectPerVoxel
      # it's used to remove voxels lesioned in few subjects
      if (!is.numeric(minSubjectPerVoxel) & is.character(minSubjectPerVoxel)) { # input is percentage
        thresholdPercent = as.numeric(gsub('%','', minSubjectPerVoxel)) / 100
      } else if (is.numeric(minSubjectPerVoxel)) { # user defined exact subject number
        thresholdPercent = minSubjectPerVoxel / length(lesions.list)
      }

      # compute average map
      avgles = antsAverageImages(lesions.list)

      # we remove voxels with too few, or too many, subjects
      mask = thresholdImage(avgles, thresholdPercent, 1 - thresholdPercent)
      # if user set minSubjectPerVoxel=0, mask is all 1, so fix it
      if (thresholdPercent == 0) mask[avgles==0] = 0
      writeavgles = TRUE

      if (showInfo) printInfo(paste(sum(mask), 'found'), type='tail') # voxels in mask
      # check mask is not empty
      if (max(mask) == 0) stop('Mask is empty. No voxels to run VLSM on.')
    }



    ##########
    # get patch information
    haslesmat = FALSE
    if (noPatch) {
      if (showInfo) printInfo('noPatch true - Patches will not be used...')
      voxmask = mask
      voxindx = 1:sum(mask)
      patchinfo.derived = 'not computed'
    } else {
      if (is.na(patchinfo[1])) {

        if (showInfo) printInfo('Computing unique patches...')
        patchinfo = getUniqueLesionPatches(lesions.list, mask=mask, showInfo = F, returnPatchMatrix = T)
        lesmat = patchinfo$patchmatrix
        haslesmat = TRUE
        patchinfo.derived = 'Computed'

      } else {

        if (showInfo) printInfo('Using predefined patch information...')
        if (!(checkImageList(list(patchinfo$patchimg, mask), showError = F)))
          stop('Patch image header is different from mask header.')

        if('patchmatrix' %in% names(patchinfo)) {
          lesmat = patchinfo$patchmatrix
          haslesmat = T
        }
        patchinfo.derived = 'Predefined'
      }

      voxindx = patchinfo$patchindx
      voxmask = patchinfo$patchimg.samples

      if (showInfo) printInfo(paste('Found',patchinfo$npatches,'patches in', patchinfo$nvoxels, 'voxels -', round(patchinfo$nvoxels/patchinfo$npatches,1), 'times more voxels'))
    }


    # create matrix of lesions
    if (showInfo & !haslesmat) printInfo('Computing lesion matrix... ', type='head')
    if (showInfo & haslesmat) printInfo('Using existing lesion matrix... ', type='head')
    if (inputtype == 'antsImageList') { # list of antsImages
      if (!haslesmat) lesmat = imageListToMatrix(lesions.list, voxmask)
    } else if (inputtype == 'antsFiles') { # vector of filenames
      if (!haslesmat) lesmat = imagesToMatrix(lesions.list, voxmask)
    }

    if (showInfo) printInfo(paste(dim(lesmat), collapse='x'), type='tail')


    ###########
    # check the matrix is binary
    # at this point discrepancies come only from unchecked filenames
    if (binaryCheck) {
      if ( !all(lesmat %in% 0:1) )
        stop('Voxels other than 0/1 detected in lesion matrix.
           To find offending image, try:
           lesions=imageFileNames2ImageList(filenames)
           checkImageList(lesions,binaryCheck=T).')
    }


    # residualize by lesion size eventually
    if ( correctByLesSize %in% c('behavior','both') ) {
      if (showInfo) printInfo('Correcting for lesion size: behavior...')
      behavior = residuals(lm(behavior ~ getLesionSize(lesions.list, showInfo)))
    }
    if ( correctByLesSize %in% c('voxel', 'both') ) {
      if (showInfo) printInfo('Correcting for lesion size: voxel...')
      lesvals = 1/sqrt(getLesionSize(lesions.list, showInfo))
      lesmat = apply(lesmat, 2, function(x) x*lesvals )
    }



    # run the analysis
    if (showInfo) printInfo(paste('Running analysis:', method,'...'), type='head')
    if (method == 'BM') {
      lsm = lsm_BM(lesmat, behavior, nperm=nperm, ...)
    } else if (method == 'BMfast') {
      lsm = lsm_BMfast(lesmat, behavior, FWERperm = (multipleComparison=='FWERperm'),
                       pThreshold = pThreshold, nperm=nperm,
                       showInfo=showInfo, ...)
    } else if (method == 'regres') {
      lsm = lsm_regres(lesmat, behavior)
    } else if (method == 'regresPerm') {
      lsm = lsm_regresPerm(lesmat, behavior)
    } else if (method == 'regresfast') {
      lsm = lsm_regresfast(lesmat, behavior,
                           FWERperm = (multipleComparison=='FWERperm'), nperm=nperm, pThreshold=pThreshold,
                           clusterPerm=(multipleComparison=='clusterPerm'), mask=mask, voxindx=voxindx, samplemask=voxmask,
                           showInfo=showInfo, ...)
    } else if (method == 'ttest') {
      lsm = lsm_ttest(lesmat, behavior, showInfo=showInfo, ...)
    } else if (method == 'welch') {
      lsm = lsm_ttest(lesmat, behavior, var.equal = FALSE, showInfo=showInfo, ...)
    } else if (method == 'chisq') {
      lsm = lsm_chisq(lesmat, behavior, runPermutations = FALSE)
    } else if (method == 'chisqPerm') {
      lsm = lsm_chisq(lesmat, behavior, runPermutations = TRUE, nperm=nperm)
    } else if (method == 'sccan') {
      lsm = lsm_sccan(lesmat, behavior, mask=mask, showInfo=showInfo, pThreshold=pThreshold, ...)
    } else if (method == 'svr') {
      lsm = lsm_svr(lesmat, behavior, showInfo=showInfo, ...)
    }
    if (showInfo) printInfo('', type='tail')


    # START POST PROCESSING THE OUTPUT
    statistic = lsm$statistic

    if ('zscore' %in% names(lsm)) {
      haszscore = T
      zscore = lsm$zscore
    } else { haszscore=F }

    if ('pvalue' %in% names(lsm)) {
      haspvalue = T
      pvalue = lsm$pvalue
    } else { haspvalue=F }


    # multiple comparison correction
    if (haspvalue) {
      if (multipleComparison %in% p.adjust.methods ) {
        if (showInfo & multipleComparison != 'none') printInfo(paste('Correcting p-values:',multipleComparison,'...'))
        pvalue.adj = p.adjust(pvalue, method = multipleComparison)
        statistic[pvalue.adj>=pThreshold] = 0
        if (haszscore) zscore[pvalue.adj>pThreshold] = 0
      } else {
        pvalue.adj = pvalue
      }
    }

    # flip sign if asked
    if (flipSign) {
      if (showInfo) printInfo('Flipping statistics negative/positive...')
      statistic = -(statistic)
      if (haszscore) zscore = -(zscore)
    }

    # a final check to make sure there are no
    # infinite, na, or nan values
    if (any(is.nan(statistic))) printInfo('WARNING: NaN values detected in statistic.', type='tail')
    if (any(is.na(statistic))) printInfo('WARNING: NA values detected in statistic.', type='tail')
    if (any(is.infinite(statistic))) printInfo('WARNING: Infinite values detected in statistic.', type='tail')
    if (haspvalue && any(is.nan(pvalue.adj))) printInfo('WARNING: NaN values detected in pvalues.', type='tail')
    if (haspvalue && any(is.na(pvalue.adj))) printInfo('WARNING: NA values detected in pvalues.', type='tail')
    if (haspvalue && any(is.infinite(pvalue.adj))) printInfo('WARNING: Infinite values detected in pvalues.', type='tail')
    if (haszscore && any(is.nan(zscore))) printInfo('WARNING: NaN values detected in zscore.', type='tail')
    if (haszscore && any(is.na(zscore))) printInfo('WARNING: NA values detected in zscore.', type='tail')
    if (haszscore && any(is.infinite(zscore))) printInfo('WARNING: Infinite values detected in zscore.', type='tail')


    # put results in images
    if (showInfo) printInfo('Preparing images...')
    vlsm.stat = makeImage(mask, voxval=statistic[voxindx])
    # optional outputs
    if (haspvalue) vlsm.pval = makeImage(mask, voxval=pvalue.adj[voxindx])
    if (haspvalue) vlsm.pval[mask==0] = 1
    if (haszscore) vlsm.zmap = makeImage(mask, voxval=zscore[voxindx])



    # call details
    if (showInfo) printInfo('Logging call details...')
    noinfo = c('lesions.list', 'mask', 'patchinfo', 'behavior', '...') # skip info on these vars
    callinfo = mget(names(formals()),sys.frame(sys.nframe()))
    callinfo = callinfo[! names(callinfo) %in% noinfo]
    mcall = as.list(match.call())[-1]
    mcall = mcall[!names(mcall) %in% noinfo]
    addindx = (! names(mcall) %in% names(callinfo))
    callinfo = c(callinfo, mcall[addindx])
    callinfo$Subjects = length(behavior)
    callinfo$patchinfo = patchinfo.derived
    callinfo = c(LesymapVersion=ver,  callinfo)
    callinfo = c(Time=as.character(toc), callinfo)




    # return results
    output = list(stat.img=vlsm.stat)
    # optional
    if (exists('vlsm.zmap')) output$zmap.img=vlsm.zmap
    if (exists('vlsm.pval')) output$pval.img=vlsm.pval

    output$mask.img=mask
    otherindx = (! names(lsm) %in% c('statistic', 'pvalue', 'zscore'))
    output = c(output, lsm[otherindx])
    if (writeavgles) output$average.img=avgles
    if (!noPatch) output$patchinfo = patchinfo
    output$callinfo = callinfo


    # save results
    if (!is.na(saveDir)) {
      if (showInfo) printInfo('Saving on disk...')
      save.lesymap(lsm = output, saveDir = saveDir, callinfo = callinfo,  ...)
    }

    tic = Sys.time()
    runtime = paste(round(as.double(difftime(tic,toc)),1), units(difftime(tic,toc)))
    output$callinfo$Runtime = runtime
    if (showInfo) printInfo(paste('Done!',runtime))

  }, split = TRUE, type = "output") # end printedOutput

  output$outputLog = outputLog

  rm(lesions.list)
  invisible(gc())
  class(output) = c('lesymap', class(output))
  return(output)

}
