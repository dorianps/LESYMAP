#' @export
#' @title Optimization of SCCAN sparseness
#'
#' @description
#' Function used to optimize SCCAN sparseness for
#' lesion to symptom mapping.
#'
#' @param lesmat lesion matrix
#' @param behavior behavior vector
#' @param mask antsImage mask
#' @param nFolds how many folds to use
#' @param sparsenessPenalty penalty term
#' @param lowerSparseness minimum searched sparseness
#' @param upperSparseness maximum searched sparseness
#' @param tol tolerance value, see \code{optimize()} in R
#' @param justValidate just check the CV of provided sparseness
#' @param cvRepetitions number of cross-validations at each sparseness
#' value. Dynamically set depending on sample size: <=30 to 6 reps,
#' <=40 to 5 reps, <=50 to 4 reps, > 50 to 3 reps.
#' @param mycoption standard SCCAN parameter
#' @param robust standard SCCAN parameter
#' @param nvecs standard SCCAN parameter
#' @param sparseness standard SCCAN parameter
#' @param cthresh standard SCCAN parameter
#' @param its standard SCCAN parameter
#' @param npermsSCCAN SCCAN permutations
#' @param smooth standard SCCAN parameter
#' @param sparseness.behav what sparsness to use for behavior
#' @param maxBased standard SCCAN parameter
#' @param directionalSCCAN (default=TRUE) switching to FALSE will
#' switch sparseness range in the positive side, 0.005 to 0.9
#' @param showInfo logical (default=TRUE) display messages
#' @param ... other arguments received from \code{\link{lesymap}}
#' or \code{\link{lsm_sccan}}.
#'
#' @return
#' List with: \cr
#' \code{minimum} - best sparseness value \cr
#' \code{objective} - minimum value of objective function \cr
#' \code{CVcorrelation} - cross-validated correlation of optimal sparness \cr
#'
#' @author Dorian Pustina
#'
#'
optimize_SCCANsparseness <- function(lesmat, behavior, mask,
                                     nFolds = 4,
                                     sparsenessPenalty=0.03,
                                     lowerSparseness=-0.9,
                                     upperSparseness=0.9,
                                     tol = 0.03,
                                     justValidate=FALSE,
                                     cvRepetitions=ifelse(length(behavior)<=30,6,
                                                   ifelse(length(behavior)<=40,5,
                                                   ifelse(length(behavior)<=50,4,
                                                                               3))),
                                     showInfo = TRUE,
                                     directionalSCCAN=TRUE,
                      mycoption=1,
                      robust=1,
                      sparseness=NA, # 0.045,
                      nvecs=1,
                      cthresh=150,
                      its=30,
                      npermsSCCAN=0,
                      smooth=0.4,
                      sparseness.behav = -0.99,
                      maxBased=FALSE,
                      ...) {

  # flip default bounds to negative eventually
  if (!directionalSCCAN) {
    upperSparseness=0.9
    lowerSparseness=0.005
    if (showInfo) printInfo(paste('\n       directionalSCCAN=FALSE, switchin default sparseness range', lowerSparseness, 'to', upperSparseness), type='middle')
  }

  # REQUIRES CARET
  # if (! 'caret' %in% rownames(installed.packages())) stop('SCCAN optimization requires the caret package. Try installing with install.packages("caret")')
  myfolds = list()
  for (i in 1:cvRepetitions) {
    myfolds[[i]] = .createFolds(behavior, nFolds) # caret::createFolds(behavior, nFolds)
  }

  cthresh = c(cthresh,0)

  # attention, assume lesmat and behavior are passed
  # already scaled


  sccan.masks=c(mask,NA)

  #' the optimization function
  #' Will run SCCAN on each training fold, compute
  #' behavior prediction on the test fold, and finally
  #' return a cross validated correlation from entire sample
  optimfun <- function(thissparse, lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
                       mycoption=mycoption, robust=robust, myfolds=myfolds, sparseness.behav=sparseness.behav,
                       maxBased=maxBased,
                       showInfo=showInfo, sparsenessPenalty=sparsenessPenalty) {

    if (showInfo) printInfo('\n', type='middle')
    if (showInfo) printInfo(paste0('        Checking sparseness ', round(thissparse,3),' '), type='head')

    CVcorr = rep(NA, length(myfolds))
    # rmse = rep(NA, length(myfolds))
    for (cvrep in 1:length(myfolds)) {
      if (showInfo) printInfo('.', type='middle') # showing dot for working cvrep

      sparseness = c( thissparse, sparseness.behav )
      behavior.predicted = rep(NA, length(behavior))
      traincorr = rep(NA, length(myfolds[[cvrep]]))
      for (i in 1:length(myfolds[[cvrep]])) {
        fold = myfolds[[cvrep]][[i]]
        trainsccan = sparseDecom2( inmatrix = list(lesmat[ -fold,],as.matrix(behavior[-fold])),
                                   inmask=sccan.masks, mycoption=mycoption,
                              robust=robust, sparseness=sparseness, nvecs=nvecs,
                              cthresh=cthresh,its=its, perms=npermsSCCAN, smooth=smooth,
                              maxBased=maxBased)

        behavior.predicted[fold] = lesmat[fold,] %*% t(trainsccan$eig1) %*% trainsccan$eig2
        traincorr[i] = abs(trainsccan$ccasummary[[1]])

        invisible(gc())
      } # and single CV loop

      CVcorr[cvrep] = abs(cor(behavior,behavior.predicted)) # need to make sure correlation is positive
      # rmse[cvrep] = sqrt(mean((behavior.predicted - behavior)^2))
    } # end repetition of CVs loop

    CVcorr = mean(CVcorr)
    # rmse = mean(rmse)

    output = 1 - ( CVcorr - (abs(thissparse)*sparsenessPenalty) )

    if (showInfo) printInfo (paste0(' CV correlation ', format(CVcorr,digits=3,nsmall=3),
                             ' (', format(mean(traincorr, na.rm=T),digits=3,nsmall=3), ')',
                             ' (cost=',format(output,digits=3,nsmall=3), ')' ), type='middle')

    return(output)
  }
  #'
  #' end of optimfun
  #'

  if (!justValidate) { # FULL OPTIMIZATION
    if (showInfo) printInfo(paste('\n       Searching for optimal sparseness:'), type='middle')
    if (showInfo) printInfo(paste('\n         lower/upper bound:\t ', lowerSparseness, '/', upperSparseness), , type='middle')
    if (showInfo) printInfo(paste('\n         cvRepetitions:\t\t ', cvRepetitions), type='middle')
    if (showInfo) printInfo(paste('\n         nFolds:\t\t ', nFolds), type='middle')
    if (showInfo) printInfo(paste('\n         sparsenessPenalty:\t ', sparsenessPenalty), type='middle')
    if (showInfo) printInfo(paste('\n         optim tolerance:\t ', tol), type='middle')

    temp = optimize(f=optimfun, lower=lowerSparseness, upper=upperSparseness, maximum = F, tol = tol,
                    lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
                    mycoption=mycoption, robust=robust, myfolds=myfolds, sparseness.behav=sparseness.behav,
                    maxBased=maxBased,
                    sparsenessPenalty=sparsenessPenalty,
                    showInfo=showInfo
                    )
    # recover true CV correlation, not just the returned function outcome
    temp$CVcorrelation.stat = (1 - temp$objective) + (abs(temp$minimum)*sparsenessPenalty)
  } else { # JUST CHECK THE CV FOR USER DEFINED SPARSENESS
    if (showInfo) printInfo('\n       Validating sparseness:', type='middle')
    if (showInfo) printInfo(paste('\n         cvRepetitions:\t\t ', cvRepetitions), type='middle')
    if (showInfo) printInfo(paste('\n         nFolds:\t\t ', nFolds), type='middle')
    objective = optimfun(thissparse=sparseness, lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
                         mycoption=mycoption, robust=robust, myfolds=myfolds, sparseness.behav=sparseness.behav,
                         maxBased=maxBased,
                         sparsenessPenalty=sparsenessPenalty,
                         showInfo=showInfo)
    temp=list()
    # recover recover true CV correlation, not just the returned function outcome
    temp$CVcorrelation.stat = (1 - objective) + (abs(sparseness)*sparsenessPenalty)
    temp$minimum = sparseness
  }


  return(temp)

}
