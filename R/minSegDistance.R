#' @title Metric displacement of binary masks
#'
#' @description
#' This function computes the metric displacement between two binary masks
#'
#' @param manual manual segmentation of class antsImage,
#' used as reference
#' @param predict other antsImage to compare to manual
#' @param get (default='all') one of 'mean', 'max', 'min', or 'all'
#' @param binarize logical (default=FALSE) whether to binarize
#' the input images
#' @param label (default=1) integer or vector of labels to binarize
#' I.e., label=c(2,4) means label 2 from manual, and 4 from predict
#' will be compared.
#'
#' @return Scalar (for 'mean', 'max', 'min') or list (for 'all').
#' Note, results are in milimeters
#'
#' @note max = Hausdorff distance
#'
#' @author Dorian Pustina
#'
#' @export
minSegDistance = function (manual, predict, get='all', binarize=F, label=1) {

  if (class(manual)!='antsImage' || class(predict)!='antsImage') stop('Images must be of class antsImage')

  if (binarize) {
    if (length(label)==2) {
      manual = thresholdImage(manual, label[1], label[1])
      predict = thresholdImage(predict, label[2], label[2])
    } else {
      manual = thresholdImage(manual, label[1], label[1])
      predict = thresholdImage(predict, label[1], label[1])
    }
  }

  pdist = predict - (predict %>% iMath('ME', 1))
  pdist = pdist %>% iMath('D')

  mbord = manual - (manual %>% iMath('ME', 1))

  if (get=='mean')
    return ( mean(pdist[mbord==1]) )
  if (get=='min')
    return ( min(pdist[mbord==1]) )
  if (get=='max')
    return ( max(pdist[mbord==1]) )
  if (get=='sum')
    return ( sum(pdist[mbord==1]) )
  if (get=='all')
    return ( list(
      mean=mean(pdist[mbord==1]),
      min=min(pdist[mbord==1]),
      max=max(pdist[mbord==1]),
      sum=sum(pdist[mbord==1])
      ) )
}
