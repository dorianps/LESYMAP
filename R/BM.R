#' @title Slow R-based Brunner-Munzel tests
#'
#' @description
#' Takes a binary matrix of voxels and a vector of behavior
#' and runs Brunner-Munzel tests on each voxel.
#' This function is not compiled and is slow.
#'
#' @param lesmat matrix of voxels
#' @param behavior vector of behavior
#' @return
#' Returned list with:
#' \itemize{
#' \item\code{statistic} statistical values
#' \item\code{dof} degrees of freedom
#' }
#'
#' @examples
#' set.seed(123)
#' lesmat = matrix(rbinom(200,1,0.5), ncol=2)
#' set.seed(123)
#' behavior = rnorm(100)
#' result = BM(lesmat, behavior)
#'
#' @author Dorian Pustina
#'
#' @export
BM <- function(lesmat, behavior) {

  statistic = dfbm = rep(NA, ncol(lesmat))

  for (i in 1:length(statistic)) {
    x=behavior[lesmat[,i]==0]
    y=behavior[lesmat[,i]==1]
    n1 = length(x)
    n2 = length(y)
    r1 = rank(x)
    r2 = rank(y)
    r = rank(c(x, y))
    m1 = mean(r[1:n1])
    m2 = mean(r[n1 + 1:n2])
    pst = (m2 - (n2 + 1)/2)/n1
    v1 = sum((r[1:n1] - r1 - m1 + (n1 + 1)/2)^2)/(n1 - 1)
    v2 = sum((r[n1 + 1:n2] - r2 - m2 + (n2 + 1)/2)^2)/(n2 - 1)
    statistic[i] = n1 * n2 * (m1 - m2)/(n1 + n2)/sqrt(n1 * v1 +  n2 * v2)
    dfbm[i] = ((n1 * v1 + n2 * v2)^2)/(((n1 * v1)^2)/(n1 - 1) + ((n2 * v2)^2)/(n2 - 1))
  }




 return(list(
    statistic=statistic,
    dof=dfbm))
}
