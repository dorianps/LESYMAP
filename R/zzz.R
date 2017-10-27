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
  gitversion = gsub('^Version: ', '', description.lines[verline[1]])
  installversion = as.character(packageVersion('LESYMAP'))
  newversion = compareVersion(gitversion,installversion)

  # display message
  if (newversion == 1) {
    packageStartupMessage(paste0('New LESYMAP version available: ', gitversion, ' (current ', installversion,')'))
    packageStartupMessage('See version history at https://git.io/vFf6g')
  }
}
