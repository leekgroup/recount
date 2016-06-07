#' This function downloads the metadata for all projects.
#'
#' Download the metadata from all the projects. This can be useful for finding
#' samples of interests across all projects.
#'
#' @param subset Either \code{sra} or \code{gtex}. Specifies which metadata
#' file to download.
#' @param verbose If \code{TRUE} it will print a message of where the file is
#' being downloaded to.
#'
#' @return A \link[S4Vectors]{DataFrame-class} object with the phenotype
#' metadata.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' metadata <- all_metadata()
#'
#'

all_metadata <- function(subset = 'sra', verbose = TRUE) {
    ## check inputs
    subset <- tolower(subset)
    stopifnot(subset %in% c('sra', 'gtex'))
    stopifnot(length(subset) == 1)
    
    ## Download file
    metafile <- paste0('metadata_clean_', subset, '.Rdata')
    destfile <- file.path(tempdir(), metafile)
    if(verbose) message(paste(Sys.time(), 'downloading the metadata to', destfile))
    download.file(paste0(
        'https://github.com/leekgroup/recount-website/blob/master/metadata/',
        metafile, '?raw=true'), destfile = destfile) 
    load(destfile)
    return(metadata_clean)
}
