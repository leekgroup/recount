#' Add additional curated metadata to a recount rse object
#'
#' This function appends sample metadata information to a
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class] from the
#' recount2 project. The sample metadata comes from curated efforts
#' independent from the original recount2 project. Currently the only
#' information comes from the recount_brain project described in more detail
#' at <http://lieberinstitute.github.io/recount-brain/>.
#'
#'
#' @param rse A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object as downloaded with [download_study]. If this argument is
#' not specified, the function will return the raw metadata table.
#' @param source A valid source name. The only supported options at this
#' moment are `recount_brain_v1` and `recount_brain_v2`.
#' @param is_tcga Set to `TRUE` only when `rse` is from TCGA.
#' Otherwise set to `FALSE` (default).
#' @param verbose If `TRUE` it will print a message of where the
#' predictions file is being downloaded to.
#'
#' @return A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object with the sample metadata columns appended to the `colData()`
#' slot.
#'
#' @details
#' For `source = "recount_brain_v1"` and
#' `source = "recount_brain_v2"`, the metadata columns are
#' described at <http://lieberinstitute.github.io/recount-brain/>.
#' Alternatively, you can explore `recount_brain_v2` interactively at
#' <https://jhubiostatistics.shinyapps.io/recount-brain/>.
#'
#' If you use the recount_brain data please cite the Razmara et al.
#' bioRxiv, 2019 <https://www.biorxiv.org/content/10.1101/618025v1>.
#' A bib file is available via citation('recount')[5].
#'
#'
#' @references
#' Razmara et al, bioRxiv, 2019.
#' <https://www.biorxiv.org/content/10.1101/618025v1>
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
#' rse_gene <- add_metadata(rse_gene_SRP009615, 'recount_brain_v2')
#'
#' ## Explore the metadata
#' colData(rse_gene)
#'
#' ## For a list of studies present in recount_brain check
#' ## http://lieberinstitute.github.io/recount-brain/.
#' ## recount_brain_v2 includes GTEx and TCGA brain samples in addition to the
#' ## recount_brain_v1 data, plus ontology information.
#'
#'
#' ## Obtain all the recount_brain_v2 metadata if you want to
#' ## explore the metadata manually
#' recount_brain_v2 <- add_metadata(source = 'recount_brain_v2')
#'

add_metadata <- function(rse, source = 'recount_brain_v2', is_tcga = FALSE,
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
    download_retry(
        url = to_use$url,
        destfile = destfile,
        mode = 'wb'
    )
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
