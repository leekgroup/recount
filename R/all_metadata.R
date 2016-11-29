#' This function downloads the metadata for all projects.
#'
#' Download the metadata from all the projects. This can be useful for finding
#' samples of interests across all projects.
#'
#' @param subset Either \code{sra},  \code{gtex} or \code{tcga}. Specifies
#' which metadata file to download.
#' @param verbose If \code{TRUE} it will print a message of where the file is
#' being downloaded to.
#'
#' @return A \link[S4Vectors]{DataFrame-class} object with the phenotype
#' metadata.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import downloader
#'
#' @examples
#' 
#' metadata <- all_metadata()
#'
#' @details Note that for \code{subset = 'gtex'}, there are more variables than
#' the ones we have for 'sra'. This information corresponds to file
#' GTEx_Data_V6_Annotations_SampleAttributesDS.txt available at
#' \url{http://www.gtexportal.org/home/datasets}. There you can find the
#' information describing these variables.
#'
#' For TCGA we acquired metadata information from 3 different sources:
#' - GDC: via a json query
#' - CGC: via json queries and a custom script to merge the tables
#' - TCGAbiolinks: we used to to parse GDC's XML files
#' For more information, check \url{https://github.com/leekgroup/recount-website/tree/master/metadata/tcga_prep}.
#'

all_metadata <- function(subset = 'sra', verbose = TRUE) {
    ## For R CMD check
    metadata_clean <- NULL
    
    ## check inputs
    subset <- tolower(subset)
    stopifnot(subset %in% c('sra', 'gtex', 'tcga'))
    stopifnot(length(subset) == 1)
    
    ## Download file
    metafile <- paste0('metadata_clean_', subset, '.Rdata')
    url <- paste0(
        'https://github.com/leekgroup/recount-website/blob/master/metadata/',
        metafile, '?raw=true')
    destfile <- file.path(tempdir(), metafile)
    
    if(verbose) message(paste(Sys.time(), 'downloading the metadata to', destfile))
    downloader::download(url, destfile = destfile, mode = 'wb') 
    load(destfile)
    return(metadata_clean)
}
