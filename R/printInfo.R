#' printInfo
#'
#' Displays in the R console the information sent
#' by LESYMAP functions.
#'
#' @param message character binary matrix (0/1) of voxels (columns)
#' and subjects (rows).
#' @param type decides how to add newlines and timestamp by choosing
#'   between the types of messages LESYMAP needs.
#'   Choose between "full", "head", "middle", "tail".
#' @param tstamp the timestamp format to display
#' @param ... other arguments that are passed by upstream functions
#'
#' @return
#' NULL
#'
#' @author Dorian Pustina
#'
#' @export

printInfo <- function(message, type = 'full', tstamp = "%H:%M:%S", ...) {

  # check the type of message is accepted
  if (! type %in% c('head','tail','middle','full')) stop('Wrong printInfo type.')

  printmsg = '' # initiate

  # head stuff
  if (type %in% c('head', 'full')) {
    printmsg = paste0(printmsg, format(Sys.time(), tstamp))
  }

  # message
  printmsg = paste(printmsg, message)

  # tail stuff
  if (type %in% c('full', 'tail')) {
    printmsg = paste0(printmsg, '\n')
  }

  cat(printmsg)
  # return()
}
