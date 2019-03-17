#' lsm_chisq
#'
#' Lesion to symptom mapping performed on a prepared matrix.
#' The behavior must be a binary vector.
#' Chi square tests are performed at each voxel. By default
#' the Yates correction is performed, use \code{correct=FALSE}
#' if you need to disable it. The behavior must
#' be a binary vector. Exact p-values can be obtained with permutation
#' based estimatins.
#'
#' @param lesmat binary matrix (0/1) of voxels (columns)
#' and subjects (rows).
#' @param behavior vector of behavioral scores (must be binary.
#' @param YatesCorrect (default=T) logical whether to use Yates correction.
#' @param runPermutations logical (default=FALSE) whether to
#' use permutation based p-value estimation.
#' @param nperm (default=2000) The number of permutations to run.
#' @param showInfo display info messagges when running the function.
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
#' set.seed(1234)
#' behavior = rbinom(100,1,0.5)
#' result = lsm_chisq(lesmat, behavior)
#' }
#'
#' @author Dorian Pustina
#'
#' @export
lsm_chisq <- function(lesmat, behavior, YatesCorrect=TRUE,
                      runPermutations = F, nperm=2000,
                      showInfo=TRUE, ...) {

  behavOn = sum(behavior == 1)
  behavOff = length(behavior) - behavOn
  lesVox = colSums(lesmat)
  lesOnBehavOn = colSums(apply(lesmat, 2, function(x) x*behavior))
  lesOnBehavOff = lesVox - lesOnBehavOn
  lesOffBehavOn = behavOn - lesOnBehavOn
  lesOffBehavOff = behavOff - lesOnBehavOff
  chimatrix = rbind( lesOnBehavOff, lesOnBehavOn, lesOffBehavOff, lesOffBehavOn)
  rm(behavOn,behavOff,lesVox,lesOnBehavOn,lesOnBehavOff,lesOffBehavOn,lesOffBehavOff)

  # estiamte runtime
  if (showInfo & runPermutations) {
    tic = Sys.time()
    temp = chisq.test(matrix(chimatrix[,1],ncol=2),
               simulate.p.value=runPermutations,
               correct=YatesCorrect, B=nperm)
    onerun = as.double(difftime(Sys.time(),tic, units = 'sec'))
    toc = tic + (ncol(chimatrix)* onerun)
    expect = paste(round(as.double(difftime(toc,tic)),1), units(difftime(toc,tic)))
    printInfo(paste('\n       Chi square permutations, expected run =', expect), type='middle')
  }


  output = apply(chimatrix, 2, function(x) {
                          temp=chisq.test(matrix(x,ncol=2),
                                          simulate.p.value=runPermutations,
                                          B=nperm)
                          return(list(
                            stat=temp$statistic,
                            pval=temp$p.value))
                        }
                )

  temp = unlist(output)
  statistic = unname( temp[seq(1,length(temp),by=2)] )
  pvalue = unname( temp[seq(2,length(temp),by=2)] )
  # zscore = qchisq(pvalue, df=1)
  # zscore = qnorm(pvalue, lower.tail=TRUE)
  # zscore[is.infinite(zscore)] = 0 # fixing infinite values for p=1

  return(list(
    statistic=statistic,
    pvalue=pvalue
    #zscore=zscore
    ))
}
