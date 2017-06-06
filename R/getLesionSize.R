#' getLesionSize
#'
#' Compute lesion sizes from a list of antsImages.
#'
#' @param lesions.list List of antsImages. For proper
#' measurement, images must be binary.
#'
#' @return vector of lesion sizes
#'
#' @author Dorian Pustina
#'
#' @export
getLesionSize <- function(lesions.list) {
  lessize = rep(NA, length(lesions.list))
  voxsize = prod(antsGetSpacing(lesions.list[[1]]))
  for (i in 1:length(lesions.list)) {
    lessize[i] = sum(lesions.list[[i]]>0) * voxsize
  }
  invisible(gc()) # release memory to help R in linux
  return(lessize)
}
