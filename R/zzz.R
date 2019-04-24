# show a message during startup
.onAttach <- function(libname, pkgname){
  packageStartupMessage("Please send bug reports and suggestions to shahrisa.arman@hotmail.com with a meaningful email subject")
}
