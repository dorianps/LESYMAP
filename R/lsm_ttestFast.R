#' @title T-tests for symptom mapping (fast)
#'
#' @description
#' Lesion to symptom mapping performed on a prepared matrix.
#' T-tests are performed using each column of the
#' matrix to split the behavioral scores in two groups. If
#' var.equal=TRUE the Welch test is performed instead.
#' This function relies on TTfast, a compiled version to run
#' on thousands of voxels (about 60 faster then regular t-tests
#' in R).
#'
#' @param lesmat binary matrix (0/1) of voxels (columns)
#' and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param var.equal logical (default=TRUE) should the
#' variance between groups considered equal (t-test) or
#' unequal (Welch test).
#' @param alternative (default='greater') Sets the expected
#' relationship between voxel value and behavior. By default
#' voxels with zero are not lesioned, and behavior is expected to
#' be higher, thus \code{alternative='greater'}. If the relationship in your
#' data is inverted, use \code{alternative='less'}, and if
#' you don't have a relationship hypothesis data,
#' use \code{alternative='two.sided'}.
#' @param checkAssumptions (default=TRUE) Check how many voxels violate
#' the t-test assumptions (heteroscadsticity and/or normality).
#' @param showInfo logical (default=TRUE), display time-stamped info messages
#' @param ... other arguments received from \code{\link{lesymap}}.
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{statistic} - vector of statistical values
#' \item\code{pvalue} - vector of pvalues
#' \item\code{zscore} - vector of zscores
#' }
#'
#' @examples{
#' set.seed(123)
#' lesmat = matrix(rbinom(200,1,0.5), ncol=2)
#' set.seed(123)
#' behavior = rnorm(100)
#' result = lsm_ttestFast(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
#' @export
lsm_ttestFast <- function(lesmat, behavior,
                      var.equal = TRUE,
                      alternative='greater',
                      checkAssumptions = TRUE,
                      showInfo = TRUE,
                      ...) {

  if (! alternative %in% c('greater','less','two.sided'))
    stop('Alternative can be one of greater/less/two.sided')

  # check assumptions only for t-test, not for welch
  if (var.equal & checkAssumptions) checkAssumptions_ttest(lesmat, behavior, showInfo=showInfo, ...)

  # run t-tests
  resultlist = TTfast(X=lesmat, Y=as.matrix(behavior),
                      computeDOF = TRUE, varEqual = var.equal)

  statistic = resultlist$statistic

  # pvalue calculation
  if ((alternative == "greater") | (alternative == "g"))
    pvalue = pt(q=statistic, df=resultlist$df, lower.tail = FALSE)
  if ((alternative == "less") | (alternative == "l"))
    pvalue = pt(q=statistic, df=resultlist$df, lower.tail = TRUE)
  if (alternative == 'two.sided')
    pvalue = pt(q=abs(statistic), df=resultlist$df, lower.tail = FALSE)*2

  # zscore calculation
  if ((alternative == "less") | (alternative == "l")) {
    zscore = qnorm(pvalue, lower.tail=TRUE)
  }
  else if ((alternative == "greater") | (alternative == "g")) {
    zscore = qnorm(pvalue, lower.tail=FALSE)
  }
  else {
    zscore = qnorm(pvalue, lower.tail=FALSE)
  }

  # convert infinite values to huge numbers
  # done in lsm_ttestFast, not fixed in old lsm_ttest
  zscore[zscore == Inf] = .Machine$double.xmax
  zscore[zscore == -Inf] = -.Machine$double.xmax

  output = list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore
  )

  return(output)
}


#' @title Check t-test assumptions at each voxel
#'
#' @description
#' Routine to test statistical assumptions are met
#' at each voxel for t-tests
#'
#' @param lesmat binary matrix (0/1) of voxels (columns)
#' and subjects (rows).
#' @param behavior vector of behavioral scores.
#' @param assumptionThreshold (default=0.05) threshold to
#' establish whether tests pass or fail. This threshold is
#' applied at each voxel individually, no corretion for
#' multiple comparison is applied.
#' @param showInfo logical (default=TRUE), display info messages
#' @param ... arguments that are passed by upstream functions
#'
#' @return
#' List of objects returned:
#' \itemize{
#' \item\code{failVarianceTest} - vector of logical values
#'    marking voxels that have different variance of behavioral
#'    scores in lesioned and non-lesioned individuals. Obtained
#'    using the \code{var.test} function.
#' \item\code{failNormalityTest} - vector of logical values
#'    marking voxels with abnormal distribution of behavioral
#'    scores either in lesioned or non-lesioned individuals.
#'    Obtained with the \code{shapiro.test} function.
#' }
#'
#' @author Dorian Pustina
#'
#' @export
#' @importFrom stats shapiro.test var.test
checkAssumptions_ttest <- function(lesmat, behavior,
                                   assumptionThreshold = 0.05,
                                   showInfo = TRUE,
                                   ...) {

  # start in new line
  if (showInfo) printInfo('', type='tail')

  # test for variance homogeneity of behavioral scores at each voxel
  if (showInfo) printInfo('    checking variance homogeneity...', type='head')

  failVarianceTest = apply(lesmat, 2, function(x) var.test(behavior[x==0], behavior[x!=0])$p.value) <= assumptionThreshold

  if (showInfo) {
    msg = paste0(sum(failVarianceTest), ' voxels failed ',
            '(', round(sum(failVarianceTest)/ncol(lesmat)*100, 0), '%)')
    printInfo(msg, type='tail')
  }


  # test for normality of distribution of the
  # behavioral score for either group (x = 0 or 1) at each voxel (lesmat column)
  if (showInfo) printInfo('    checking distribution normality...', type='head')

  failNormalityTest = apply(lesmat, 2, function(x)
    (
      shapiro.test(behavior[x==0])$p.value <= assumptionThreshold |
      shapiro.test(behavior[x!=0])$p.value <= assumptionThreshold
    )
  )

  if (showInfo) {
    msg = paste0( sum(failNormalityTest), ' voxels failed ',
      '(', round(sum(failNormalityTest)/ncol(lesmat)*100, 0), '%)')
    printInfo(msg, type='middle')
  }

  # output = list()
  # output$failVarianceTest = failVarianceTest
  # output$failNormalityTest = failNormalityTest

  return()

}
