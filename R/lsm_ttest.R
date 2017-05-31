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
