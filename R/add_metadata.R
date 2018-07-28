#' Add additional curated metadata to a recount rse object
#'
#' This function appends sample metadata information to a 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} from the
#' recount2 project. The sample metadata comes from curated efforts
#' independent from the original recount2 project. Currently the only
#' information comes from the recount_brain project described in more detail
#' at \url{http://lieberinstitute.github.io/recount-brain/}.
#'
#'
#' @param rse A \link[SummarizedExperiment]{RangedSummarizedExperiment-class} 
#' object as downloaded with \link{download_study}. If this argument is
#' not specified, the function will return the raw metadata table.
#' @param source A valid source name. The only supported options at this
#' moment are \code{recount_brain_v1} and \code{recount_brain_v2}.
#' @param is_tcga Set to \code{TRUE} only when \code{rse} is from TCGA.
#' Otherwise set to \code{FALSE} (default).
#' @param verbose If \code{TRUE} it will print a message of where the 
#' predictions file is being downloaded to.
#'
#' @return A \link[SummarizedExperiment]{RangedSummarizedExperiment-class} 
#' object with the sample metadata columns appended to the \code{colData()}
#' slot.
#' For \code{source = "recount_brain_v1"}, the metadata columns are
#' described at \url{http://lieberinstitute.github.io/recount-brain/}.
#' For \code{source = "recount_brain_v2"}, the metadata columns are
#' described at
#'  \url{http://lieberinstitute.github.io/recount-brain/cross_studies_metadata/cross_studies_metadata.html}.
#'
#' @details If you use the recount_brain data please cite the Razmara et al
#' bioRxiv pre-print available at 
#' (TODO update URL once it's available). See citation
#' details with citation('recount').
#'
#' @references
#' Razmara et al, in prep, 2018.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import downloader
#' @import SummarizedExperiment
#'
#' @examples
#'
#' ## Add the sample metadata to an example rse_gene object
#' rse_gene <- add_metadata(rse_gene_SRP009615, 'recount_brain_v1')
#'
#' ## Explore the metadata
#' colData(rse_gene)
#'
#' ## For a list of studies present in recount_brain_v1 check
#' ## http://lieberinstitute.github.io/recount-brain/. Note that it only
#' ## includes studies from SRA, so no TCGA or GTEx (those have great
#' ## sample metadata already available).
#' ## recount_brain_v2 includes GTEx and TCGA brain samples in addition to the
#' ## recount_brain_v1 data.
#'
#' \dontrun{
#' ## Example project that is present in recount_brain_v2.
#'
#' ## Download and load the data
#' download_study('ERP001304')
#' load(file.path('ERP001304', 'rse_gene.Rdata'))
#'
#' ## Add the sample metadata from recount_brain_v2
#' rse_gene <- add_metadata(rse_gene, source = 'recount_brain_v2')
#'
#' ## Check the metadata
#' colData(rse_gene)
#' }
#'
#' ## Obtain all the recount_brain_v2 metadata
#' recount_brain_v2 <- add_metadata(source = 'recount_brain_v2')
#'

add_metadata <- function(rse, source = 'recount_brain_v1', is_tcga = FALSE,
    verbose = TRUE) {
        
    stopifnot(length(source) == 1)
        
    ## For a NOTE in R CMD check
    valid_sources <- data.frame(
        name = c('recount_brain_v1', 'recount_brain_v2'),
        url = c(
            'https://github.com/LieberInstitute/recount-brain/blob/master/merged_metadata/recount_brain_v1.Rdata?raw=true', 'https://github.com/LieberInstitute/recount-brain/blob/master/cross_studies_metadata/recount_brain_v2.Rdata?raw=true'),
        object = c('recount_brain', 'recount_brain'),
        sample_id = c('run_s', 'run_s'),
        stringsAsFactors = FALSE
    )
    
    stopifnot(tolower(source) %in% tolower(valid_sources$name))
    
    to_use <- valid_sources[tolower(valid_sources$name) == tolower(source), ]
    
    destfile <- file.path(tempdir(), paste0(to_use$name, '.Rdata'))
    
    
    if(verbose)  message(paste(Sys.time(), 'downloading the', to_use$object, 'metadata to', destfile))
    downloader::download(to_use$url, destfile = destfile, mode = 'wb')
    load_meta <- function() {
        load(destfile, verbose = verbose)
        get(to_use$object)
    }
    new_meta <- load_meta()
    
    if(missing(rse)) return(new_meta)
        
    if(is_tcga) {
        map <- match(colData(rse)$gdc_file_id, new_meta[, to_use$sample_id])
    } else {
        map <- match(colData(rse)$run, new_meta[, to_use$sample_id])
    }
    if(verbose) {
        message(paste(Sys.time(), 'found', sum(!is.na(map)), 'out of', length(map), 'samples in the', to_use$object, 'metadata'))
    }
    
    ## Make a dummy table with the new metadata to be added
    dummy <- as.data.frame(matrix(NA, nrow = ncol(rse),
        ncol = ncol(new_meta) - 1))
    cols_to_drop <- which(colnames(new_meta) == to_use$sample_id)
    colnames(dummy) <- colnames(new_meta)[- cols_to_drop]
    
    ## In case new data is present
    if(any(!is.na(map))){
        dummy[!is.na(map), ] <- new_meta[map[!is.na(map)], - cols_to_drop]
    }    
    rownames(dummy) <- NULL
   
    ## Merge new metadata and return the rse
    colData(rse) <- cbind(colData(rse), dummy)
    return(rse)
}
