#' @title T-tests for symptom mapping (slow)
#'
#' @description
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
  if (var.equal & checkAssumptions) checkAssumptions_ttest(lesmat, behavior, showInfo=showInfo, ...)

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

