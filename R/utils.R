#' Attempt to load the namespace of a package and install it if it's missing
#'
#' This function uses requireNamespace to try to load a package. But if it's
#' misssing it will then install it via Bioconductor.
#'
#' @param pkg A single character vector with the name of the package.
#' @param quietly Whether to run requireNamespace and BiocManager::install
#' quietly or not.
#'
#' @details Taken from the regionReport package
#'
#' @author Leonardo Collado-Torres
#'
.load_install <- function(pkg, quietly = TRUE) {
    attemptName <- requireNamespace(pkg, quietly = quietly)
    if(!attemptName) {
        attemptInstall <- tryCatch(BiocManager::install(pkg,
            suppressUpdates = quietly),
            warning = function(w) 'failed')
        if(attemptInstall == 'failed') stop(paste('Failed to install', pkg))
        attemptName <- requireNamespace(pkg, quietly = quietly)
    }
    if(attemptName) {
        if(quietly) {
            suppressPackageStartupMessages(library(package = pkg,
                character.only = TRUE))
        } else {
            library(package = pkg, character.only = TRUE)
        }
    }
    return(invisible(NULL))
}
