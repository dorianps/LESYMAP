#' lsm_ttest
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' T-tests are performed using each column of the
#' matrix to split the behavioral scores in two groups. If
#' var.equal=TRUE the Welch test is performed instead.
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
#' @param checkAssumptions Check whether t-test assumptions are met
#'   for every voxel
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
#' result = lsm_ttest(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
#' @export
lsm_ttest <- function(lesmat, behavior,
                      var.equal = T,
                      alternative='greater',
                      checkAssumptions = TRUE,
                      showInfo = TRUE,
                      ...) {

  # check assumptions only for t-test, not for welch
  if (var.equal) checkAssumptions_ttest(lesmat, behavior, showInfo=showInfo, ...)

  # run t-tests
  resultlist = apply(lesmat, 2, function(x) {
                          temp=t.test(behavior[x==0], behavior[x!=0], var.equal=var.equal, alternative=alternative)
                          return(list(
                            stat=temp$statistic,
                            pval=temp$p.value))
                        }
                )

  temp = unlist(resultlist)
  statistic = temp[seq(1,length(temp),by=2)]
  pvalue = temp[seq(2,length(temp),by=2)]

  if ((alternative == "less") | (alternative == "l")) {
    zscore = qnorm(pvalue, lower.tail=TRUE)
  }
  else if ((alternative == "greater") | (alternative == "g")) {
    zscore = qnorm(pvalue, lower.tail=FALSE)
  }
  else {
    zscore = qnorm(pvalue, lower.tail=FALSE)
  }

  # divide pvalues in half to get single-tailed threshold (careful on the assumption)

  output = list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore
  )

  return(output)
}


#' checkAssumptions_ttest
#'
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

  # test for variance homogeneity of behavioral scores at each voxel
  if (showInfo) printInfo('    checking variance homogeneity...', nlStart=TRUE, nlEnd=FALSE)

  failVarianceTest = apply(lesmat, 2, function(x) var.test(behavior[x==0], behavior[x!=0])$p.value) <= assumptionThreshold

  if (showInfo) {
    msg = paste0(sum(failVarianceTest), ' voxels failed ',
            '(', round(sum(failVarianceTest)/ncol(lesmat)*100, 0), '%)')
    printInfo(msg, nlStart=FALSE, nlEnd=TRUE, showTime=FALSE)
  }


  # test for normality of distribution of the
  # behavioral score for either group (x = 0 or 1) at each voxel (lesmat column)
  if (showInfo) printInfo('    checking distribution normality...', nlStart=FALSE, nlEnd=FALSE)

  failNormalityTest = apply(lesmat, 2, function(x)
    (
      shapiro.test(behavior[x==0])$p.value <= assumptionThreshold |
      shapiro.test(behavior[x!=0])$p.value <= assumptionThreshold
    )
  )

  if (showInfo) {
    msg = paste0( sum(failNormalityTest), ' voxels failed ',
      '(', round(sum(failNormalityTest)/ncol(lesmat)*100, 0), '%)')
    printInfo(msg, nlStart=FALSE, nlEnd=FALSE, showTime=FALSE)
  }

  # output = list()
  # output$failVarianceTest = failVarianceTest
  # output$failNormalityTest = failNormalityTest

  return()

}


