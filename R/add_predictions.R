#' Add predicted phenotypes to a recount rse object
#'
#' Shannon Ellis et al (2017) predicted phenotypes based on expression data for
#' the samples in the recount2 project. Using this function you can add the
#' predictions to a
#' [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class] object
#' to the `colData()` slot.
#'
#' @param rse A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object as downloaded with [download_study]. If this argument is
#' not specified, the function will return the full predictions table.
#' @param is_tcga Set to `TRUE` only when `rse` is from TCGA.
#' Otherwise set to `FALSE` (default).
#' @param version The version number for the predicted phenotypes data. It has
#' to match one of the available numbers at
#' <https://github.com/leekgroup/recount-website/blob/master/predictions/>.
#' Feel free to check if there is a newer version than the default. The version
#' used is printed as part of the file name.
#' @param verbose If `TRUE` it will print a message of where the
#' predictions file is being downloaded to.
#'
#' @return A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object with the prediction columns appended to the `colData()` slot.
#' The predicted phenotypes are:
#' \describe{
#' \item{sex }{ male or female,}
#' \item{samplesource }{ cell_line or tissue,}
#' \item{tissue }{ tissue predicted based off of 30 tissues in GTEx,}
#' \item{sequencingstrategy }{ single or paired end sequencing.}
#' }
#' For each of the predicted phenotypes there are several columns as described
#' next:
#' \describe{
#' \item{reported_phenotype }{ `NA` when not available,}
#' \item{predicted_phenotype }{ `NA` when we did not predict, "Unassigned"
#' when prediction was ambiguous,}
#' \item{accuracy_phenotype }{ accuracy is assigned per dataset based on
#' comparison to samples for which we had reported phenotype information so
#' there are three distinct values per predictor (GTEx, TCGA, SRA) across all
#' studies.}
#' }
#'
#' @details If you use these predicted phenotypes please cite the Ellis et al
#' bioRxiv pre-print available at
#' https://www.biorxiv.org/content/early/2017/06/03/145656. See citation
#' details with citation('recount').
#'
#' @references
#' Ellis et al, bioRxiv, 2017.
#' https://www.biorxiv.org/content/early/2017/06/03/145656
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import downloader
#' @import SummarizedExperiment
#'
#' @examples
#'
#' ## Add the predictions to an example rse_gene object
#' rse_gene <- add_predictions(rse_gene_SRP009615)
#'
#' ## Explore the predictions
#' colData(rse_gene)
#'
#' ## Download all the latest predictions
#' PredictedPhenotypes <- add_predictions()
#'

add_predictions <- function(rse, is_tcga = FALSE, version = 'latest',
    verbose = TRUE) {

    ## For a NOTE in R CMD check
    PredictedPhenotypes <- NULL
    if(version == 'latest') version <- tryCatch(suppressWarnings(readLines(
        'https://raw.githubusercontent.com/leekgroup/recount-website/master/predictions/latestVersion.txt'))[1]
        , error = function(e) {
            print(e)
            v <- '0.0.05'
            message(paste(Sys.time(), 'Failed to check the latest version, using version', v))
            return(v) })

    ## Download file
    predfile <- paste0('PredictedPhenotypes_v', version, '.rda')
    url <- paste0(
        'https://github.com/leekgroup/recount-website/blob/master/predictions/',
        predfile, '?raw=true')
    destfile <- file.path(tempdir(), predfile)

    if(verbose) message(paste(Sys.time(), 'downloading the predictions to', destfile))
    download_retry(
        url = url,
        destfile = destfile,
        mode = 'wb'
    )
    load(destfile, verbose = verbose)
    if(missing(rse)) return(PredictedPhenotypes)

    if(is_tcga) {
        map <- match(colData(rse)$gdc_file_id, PredictedPhenotypes$sample_id)
    } else {
        map <- match(colData(rse)$run, PredictedPhenotypes$sample_id)
    }
    colData(rse) <- cbind(colData(rse), PredictedPhenotypes[map,
        !colnames(PredictedPhenotypes) %in% c('sample_id', 'dataset')])
    return(rse)
}
