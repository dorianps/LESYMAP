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
#' @param nfolds how many folds to use
#' @param sparsenessPenalty penalty term
#' @param lower minimum searched sparseness
#' @param upper maximum searched sparseness
#' @param tol tolerance value, see optimize() in R
#' @param justValidate just check the CV of provided sparseness
#' @param cvRepetitions number of cross-validations
#' @param mycoption standard SCCAN parameter
#' @param robust standard SCCAN parameter
#' @param nvecs standard SCCAN parameter
#' @param sparseness standard SCCAN parameter
#' @param cthresh standard SCCAN parameter
#' @param its standard SCCAN parameter
#' @param npermsSCCAN SCCAN permutations
#' @param smooth standard SCCAN parameter
#' @param sparseness.behav what sparsness to use for behavior
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
optimize_SCCANsparseness <- function(lesmat, behavior, mask, nfolds = 4, sparsenessPenalty=0.03,
                                     lower=0.005, upper=0.9, tol = 0.03, justValidate=F, cvRepetitions=3,
                      mycoption=1,
                      robust=1,
                      sparseness=NA, # 0.045,
                      nvecs=1,
                      cthresh=150,
                      its=30,
                      npermsSCCAN=0,
                      smooth=0.4,
                      sparseness.behav = -0.99,
                      ...) {

  # REQUIRES CARET
  if (! 'caret' %in% rownames(installed.packages())) stop('SCCAN optimization requires the caret package. Try installing with install.packages("caret")')
  myfolds = list()
  for (i in 1:cvRepetitions) myfolds[[i]] = caret::createFolds(behavior, nfolds)

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
                       showInfo=T, tstamp="%H:%M:%S", sparsenessPenalty=sparsenessPenalty) {

    if (showInfo) cat(paste0('\n', format(Sys.time(), tstamp), '            Checking sparseness ', round(thissparse,3),' ... '))

    CVcorr = rep(NA, length(myfolds))
    for (cvrep in 1:length(myfolds)) {

      sparseness = c( thissparse, sparseness.behav )
      behavior.predicted = rep(NA, length(behavior))
      traincorr = rep(NA, length(myfolds[[cvrep]]))
      for (i in 1:length(myfolds[[cvrep]])) {
        fold = myfolds[[cvrep]][[i]]
        trainsccan = sparseDecom2( inmatrix = list(lesmat[ -fold,],as.matrix(behavior[-fold])),
                                   inmask=sccan.masks, mycoption=mycoption,
                              robust=robust, sparseness=sparseness, nvecs=nvecs,
                              cthresh=cthresh,its=its, perms=npermsSCCAN, smooth=smooth )

        behavior.predicted[fold] = lesmat[fold,] %*% t(trainsccan$eig1) %*% trainsccan$eig2
        traincorr[i] = trainsccan$ccasummary[[1]]

        invisible(gc())
      } # and single CV loop

      CVcorr[cvrep] = cor(behavior,behavior.predicted)
    } # end repetition of CVs loop

    CVcorr = mean(CVcorr)

    output = 1 - ( CVcorr - (thissparse*sparsenessPenalty) )

    if (showInfo) cat(paste0('cross-validated correlation ', format(CVcorr,digits=3,nsmall=3),
                             ' (', format(mean(traincorr, na.rm=T),digits=3,nsmall=3), ')',
                             ' (cost=',format(output,digits=3,nsmall=3),')' ))

    return(output)
  }
  #'
  #' end of optimfun
  #'

  if (!justValidate) { # FULL OPTIMIZATION
    temp = optimize(f=optimfun, lower=lower, upper=upper, maximum = F, tol = tol,
                    lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
                    mycoption=mycoption, robust=robust, myfolds=myfolds, sparseness.behav=sparseness.behav,
                    sparsenessPenalty=sparsenessPenalty,
                    showInfo=T, tstamp="%H:%M:%S"
                    )
    # recover true CV correlation, not just the returned function outcome
    # temp$CVcorrelation.stat = 1 - temp$objective
    temp$CVcorrelation.stat = (1 - temp$objective) + (temp$minimum*sparsenessPenalty)
  } else { # JUST CHECK THE CV FOR USER DEFINED SPARSENESS
    objective = optimfun(thissparse=sparseness, lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
                         mycoption=mycoption, robust=robust, myfolds=myfolds, sparseness.behav=sparseness.behav,
                         sparsenessPenalty=sparsenessPenalty,
                         showInfo=T, tstamp="%H:%M:%S")
    temp=list()
    # recover recover true CV correlation, not just the returned function outcome
    # temp$CVcorrelation.stat = 1 - objective
    temp$CVcorrelation.stat = (1 - objective) + (sparseness*sparsenessPenalty)
    temp$minimum = sparseness
  }

#   temp = nlm(f=optimfun, p=0.1, fscale=-1, ndigit=3, gradtol=0.01, stepmax=0.2, steptol=0.01, iterlim=20,
#                   lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
#                   mycoption=mycoption, robust=robust, myfolds=myfolds, showInfo=T, tstamp="%H:%M:%S"
#   )

#   clist=list(eval.max=20, abs.tol=0.01, rel.tol=0.005, x.tol=0.005, step.min=0.01, step.max=0.2)
#   temp = nlminb(start=0.04, objective=optimfun, lower=0.01, upper=0.8, control=clist,
#              lesmat=lesmat, behavior=behavior, sccan.masks=sccan.masks, cthresh=cthresh,
#              mycoption=mycoption, robust=robust, myfolds=myfolds, showInfo=T, tstamp="%H:%M:%S"
#   )

  return(temp)

}
