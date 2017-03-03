#these are simple function that setup onAttach and onLoad calls.

.onAttach<- function(libname, pkgname) {
    packageStartupMessage("Welcome to ecomix, let's get mixing. To finite mixture models to ecological data see see ?ecomix")
}

.onLoad <-  function(libname, pkgname){
    # Generic DLL loader
    dll.path <- file.path( libname, pkgname, 'libs')
    if (nzchar( subarch <- .Platform$r_arch))
      dll.path <- file.path( dll.path, subarch)
      this.ext <- paste(sub( '.', '[.]', .Platform$dynlib.ext, fixed=TRUE), '$', sep='')
      dlls <- dir(dll.path, pattern=this.ext, full.names=FALSE)
      names( dlls) <- dlls
    if (length(dlls))
      lapply(dlls, function( x) library.dynam( sub( this.ext, '', x), package=pkgname, lib.loc=libname))
  }
