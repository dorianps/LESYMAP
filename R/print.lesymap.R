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
  cat("Lesymap Summary:\n")
  cat("  Subjects            :", x$callinfo$Subjects, "\n")
  cat("  Method              :", x$callinfo$method, "\n")
  cat("  Voxels in mask      :", sum(x$mask.img), "\n")
  cat("  Method              :", x$callinfo$method, "\n")
  cat("  Multiple comparison :", x$callinfo$multipleComparison, "\n")
  cat("  P-threshold         :", x$callinfo$pThreshold, "\n")
  cat("  Statistic range     :", paste(round(range(x$stat.img),1), collapse=" "), "\n")
  cat("  Significant voxels  :", sum(x$stat.img!=0), "\n")
  cat("  Runtime             :", x$callinfo$Runtime, "\n")
  cat("  Version             :", x$callinfo$LesymapVersion, "\n")
  cat("  ... save the output for more info ...\n")
  cat("\n")
}
