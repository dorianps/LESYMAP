#' @title Simulation of behavior scores from lesion maps
#'
#' @description
#' Function simulate behavioral scores based
#' on the lesion load of specific brain areas.
#' Used to run simulation studies.
#'
#' @param lesions.list list of lesions (antsImages)
#'  or vector of filenames.
#' @param parcellation mask or parcellation image. If a
#' parcellation is passed, lesion load will be computed
#'  for each different label (value) in the image.
#'  Zero and non-affected labels are not returned
#'  by default. The parcellation input can be an
#'  antsImage or a character vector pointing to a
#'  file.
#' @param label (default=NA) if a parcellation scheme
#' id being used, you can select which labels
#'  to simulate behaviors for (i.e.,
#'  c(101,43) to simulate behavior for
#'  labels with value 101 and 43 only). If not set
#'  a simulation will be returned from each parcel.
#' @param mask mask to restrict the count of lesioned voxels.
#'  It is not recommended to use a mask, because
#'  lesions should affect behavior as they are,
#'  without the user restricting the lesions to
#'  masks defined in post-processing.
#' @param errorWeight (default=0.5) the amount of error to
#' be added, i.e., 0.5 means half of the simulation
#'  will be error, the other half signal
#' @param binaryCheck (default=FALSE) check to make sure all lesions
#' are binary
#' @param exponent power exponent to elevate behavior in order
#'  to increase non-linearity relationship with
#'  lesion load. 1 is default, and 3 is what
#'  \href{https://www.ncbi.nlm.nih.gov/pubmed/24339811}{Wang (2013)}
#'  reported as lesion load relationship with
#'  behavior.
#'
#' @return List of objectas returned:
#' \itemize{
#'  \item\code{behavload} - a matrix of simulated behavioral scores.
#'  Each column shows simulation for a single parcel. Column names
#'  indicate the label number in the parcellation file.
#'  \item\code{lesload} - same as behavload, but indicates lesions
#'  loads of the simulated regions.
#'  \item\code{lesbehavCorrelation} - vector of Pearson correlations
#'  between lesion load and simulated scores.
#'  \item\code{LesvolBehavCorrelation} - vector of Pearson
#'  correlations between lesion size and simulated scores.
#'  }
#'
#' @examples{
#'  \dontrun{
#'   lesydata = file.path(find.package('LESYMAP'),'extdata')
#'   parcellation = antsImageRead(
#'   Sys.glob(file.path(lesydata, 'template', 'Parcellation_403areas.nii.gz')))
#'   filenames = Sys.glob(file.path(lesydata, 'lesions', '*.nii.gz'))
#'   lesions = imageFileNames2ImageList(filenames)
#'   simBehavior = simulateBehavior(lesions.list = lesions, parcellation =
#'   parcellation, label = c(101,43))
#'  }
#' }
#'
#' @author Dorian Pustina
#'
#' @export
simulateBehavior <- function(lesions.list, parcellation, label=NA,
                             mask=NA, errorWeight=0.5, binaryCheck=FALSE,
                             exponent=1) {

  if (!(errorWeight >= 0 & errorWeight <= 1)) stop('errorWeight must be between 0 and 1')

  lesload = getLesionLoad(lesions.list, parcellation,
                          binaryCheck=binaryCheck, label=label, mask=mask)


  errormat = matrix(
    runif(length(lesload),min = 0, max = 1),
    ncol=ncol(lesload)
  )


  behavWeight=1-errorWeight
  lesload = -lesload
  behavload = (errormat*errorWeight) + (lesload*behavWeight)

  behavload = behavload^exponent

  corrLesSize = cor(behavload, getLesionSize(lesions.list))

  return(list(
    behavLoad = behavload,
    lesLoad = -lesload,
    lesBehavCorrelation = diag(cor(lesload,behavload)),
    LesvolBehavCorrelation = corrLesSize
    ))

}
