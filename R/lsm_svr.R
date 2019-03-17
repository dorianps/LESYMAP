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
#' @param ... other arguments received from \code{\link{lesymap}}.
#'
#' @return
#' List of objects returned:
#' \itemize{
#'   \item\code{statistic} - vector of statistical values
#'   \item\code{pvalue} - vector of pvalues
#' }
#'
#' @export
#' @author Daniel Wiesen, Dorian Pustina

lsm_svr <- function(lesmat, behavior,
                    SVR.nperm = 10000,
                    SVR.type = 'eps-regression',
                    SVR.kernel = 'radial',
                    SVR.gamma = 5,
                    SVR.cost = 30,
                    SVR.epsilon = 0.1,
                    showInfo = TRUE,
                    ...) {


  if (! 'e1071' %in% rownames(installed.packages())) stop('SVR-LSM requires e1071 package. Try installing with install.packages("e1071")')


  # scale and center data
  behavior = scale(behavior, scale=T, center=T)
  lesmat = scale(lesmat, scale=T, center=T)

  # compute SVR
  if (showInfo) {
    printInfo(paste('\n       Calling SVR with:'), type='middle')
    printInfo(paste('\n            SVM type =', SVR.type), type='middle')
    printInfo(paste('\n            Kernel =', SVR.kernel), type='middle')
    printInfo(paste('\n            Gamma =', SVR.gamma), type='middle')
    printInfo(paste('\n            Cost =', SVR.cost), type='middle')
    printInfo(paste('\n            Epsilon =', SVR.epsilon), type='middle')
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
    printInfo(paste('\n       SVR', SVR.nperm, 'permutations, expected run =', expect), type='middle')
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

}

