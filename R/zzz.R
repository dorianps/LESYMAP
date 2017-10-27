.onAttach = function(libname = find.package("LESYMAP"),
                     pkgname = "LESYMAP")
{

  # set URL
  description = 'https://raw.githubusercontent.com/dorianps/LESYMAP/master/DESCRIPTION'

  # Try to establish a connection
  suppressWarnings( conn <- try( url(description) , silent=TRUE ) )

  # If connection, try to parse values, otherwise return NULL
  if ( all( class(conn) != "try-error") ) {
    suppressWarnings( description.lines <- try( readLines(conn) , silent=TRUE ) )
    close(conn)
  } else {
    return(NULL)
  }

  # Extract version info
  verline = grep('^Version: ', description.lines)
  if (length(verline) == 0) return(NULL) # did not find any version line
  gitversion = gsub('^Version: ', '', description.lines[verline[1]])
  installversion = as.character(packageVersion('LESYMAP'))
  newversion = utils::compareVersion(gitversion,installversion)

  # display message
  if (newversion == 1) {
    packageStartupMessage(paste0('New LESYMAP available: ', gitversion, ' (installed ', installversion,')'))
    packageStartupMessage('See changes at https://git.io/vFf6g')
  }
}
