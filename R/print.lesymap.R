#' print.lesymap
#'
#' Funciton to display some meaningful summary
#' when a lesymap output is called in command line.
#'
#' @param x the output from a lesymap() call.
#' @param ... useless for compatibility with default print.
#'
#' @author Dorian Pustina
#'
#' @export
print.lesymap <- function(x, ...) {
  printInfo("Lesymap Summary:", type='tail')
  printInfo(paste0("  Version             :", x$callinfo$LesymapVersion), type='tail')
  printInfo(paste0("  Subjects            :", x$callinfo$Subjects), type='tail')
  printInfo(paste0("  Method              :", x$callinfo$method), type='tail')
  printInfo(paste0("  Voxels in mask      :", sum(x$mask.img)), type='tail')
  printInfo(paste0("  Multiple comparison :", x$callinfo$multipleComparison), type='tail')
  printInfo(paste0("  P-threshold         :", x$callinfo$pThreshold), type='tail')
  printInfo(paste0("  Statistic range     :", paste(round(range(x$stat.img),1), collapse=" ")), type='tail')
  printInfo(paste0("  Significant voxels  :", sum(x$stat.img!=0)), type='tail')
  printInfo(paste0("  Runtime             :", x$callinfo$Runtime), type='tail')
  printInfo("  ... save the output for more info ...", type='tail')
  printInfo("", type='tail')
}
