# these are simple function that setup onAttach and onLoad calls.

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Welcome to ecomix, let's get mixing. To run finite mixture models on ecological data see ?species_mix or ?region_mix")
}

.onLoad <- function(libname, pkgname) {
    # Generic DLL loader
    dll.path <- file.path(libname, pkgname, "libs")

    if (nzchar(subarch <- .Platform$r_arch)) {
        dll.path <- file.path(dll.path, subarch)
    }

    this.ext <- paste(sub(".", "[.]", .Platform$dynlib.ext, fixed = TRUE), "$", sep = "")

    dlls <- dir(dll.path, pattern = this.ext, full.names = FALSE)

    names(dlls) <- dlls

    if (length(dlls)){
              lapply(dlls, function(x) library.dynam(sub(this.ext, "", x), package = pkgname, lib.loc = libname))
      }
}

# to use magrittr shortcut
utils::globalVariables(".")

# magrittr like functions to return something else if condition is not met
return_if_not <- function(x, test, y) {
  if (test) y else x
}

# bind y to x if y not in z
bind_if_not_in <- function(x, y, z, out=base::get(z, parent.frame())) {
  x %>%
    return_if_not(
      y %>%
        magrittr::extract2(z) %>%
        base::is.null(.),
      x %>%
        magrittr::inset2(z, out)
    )
}


