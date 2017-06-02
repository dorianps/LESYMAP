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
#' you don't know your data, use \code{alternative='two.sided'}.
#' @param ... other arguments received from \code{\link{lesymap}}.
#'
#' @return List with vectors of statistic, pvalue, and zscore.
#'
#' @author Dorian Pustina
#'
#' @export
lsm_ttest <- function(lesmat, behavior, var.equal = T, alternative='greater', ...) {

  output = apply(lesmat, 2, function(x) {
                          temp=t.test(behavior[x==0], behavior[x!=0], var.equal=var.equal, alternative=alternative)
                          return(list(
                            stat=temp$statistic,
                            pval=temp$p.value))
                        }
                )

  temp = unlist(output)
  statistic = temp[seq(1,length(temp),by=2)]
  pvalue = temp[seq(2,length(temp),by=2)]
  zscore = qnorm(pvalue)

  # divide pvalues in half to get single-tailed threshold (careful on the assumption)

  return(list(
    statistic=statistic,
    pvalue=pvalue,
    zscore=zscore
    ))
}
