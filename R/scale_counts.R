#' Scale the raw counts provided by the recount project
#'
#' In preparation for a differential expression analysis, you will have to
#' choose how to scale the raw counts provided by the recount project. Note that
#' the raw counts are the sum of the base level coverage so you have to take
#' into account the read length or simply the total coverage for the given
#' sample (default option). You might want to do some further scaling to take
#' into account the gene or exon lengths.
#'
#' @param rse A \link[SummarizedExperiment]{RangedSummarizedExperiment-class} 
#' object as downloaded with \link{download_study}.
#' @param by Either \code{auc} or \code{mapped_reads}. If set to \code{auc} it 
#' will scale the counts by the total coverage of the sample. That is, the area
#' under the curve (AUC) of the coverage. If set to \code{mapped_reads} it will
#' scale the counts by the number of mapped reads, whether the library was
#' paired-end or not, and the desired read length (\code{L}).
#' @param targetSize The target library size in number of single end reads.
#' @param L The target read length. Only used when \code{by = 'mapped_reads'}
#' since it cancels out in the calculation when using \code{by = 'auc'}.
#' @param factor_only Whether to only return the numeric scaling factor or
#' to return a \link[SummarizedExperiment]{RangedSummarizedExperiment-class} 
#' object with the counts scaled. If set to \code{TRUE}, you have to multiply 
#' the sample counts by this scaling factor.
#' @param round Whether to round the counts to integers or not.
#'
#' @return If \code{factor_only = TRUE} it returns a numeric vector with the
#' scaling factor for each sample. If \code{factor_only = FALSE} it returns a
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object with
#' the counts already scaled.
#'
#' @details Rail-RNA \url{http://rail.bio} uses soft clipping when aligning
#' which is why we recommed using \code{by = 'auc'}.
#'
#' If the reads are from a paired-end library, then the \code{avg_read_length}
#' is the average fragment length. This is taken into account when using
#' \code{by = 'mapped_reads'}.
#'
#' @seealso \link{download_study}
#'
#' @importFrom methods is
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' ## Load an example rse_gene object
#' rse_gene <- rse_gene_SRP009615
#'
#' ## Scale counts
#' rse <- scale_counts(rse_gene)
#'
#' ## Find the project used as an example
#' project_info <- abstract_search('GSE32465')
#'
#' ## See some summary information for this project
#' project_info
#'
#' ## Use the following code to re-download this file
#' \dontrun{
#' ## Download
#' download_study(project_info$project)
#' 
#' ## Load file
#' load(file.path(project_info$project, 'rse_gene.Rdata'))
#' identical(rse_gene, rse_gene_SRP009615)
#' }


scale_counts <- function(rse, by = 'auc', targetSize = 4e7, L = 100,
    factor_only = FALSE, round = TRUE) {
        
    ## For R CMD check
    assay <- colData <- 'assay<-' <- NULL
    
    ## Load required packages
    .load_install('SummarizedExperiment')
    .load_install('GenomicRanges')
    
    ## Check input
    stopifnot(is(rse, 'RangedSummarizedExperiment'))
    stopifnot(length(targetSize) == 1)
    stopifnot(is.numeric(targetSize) | is.integer(targetSize))
    stopifnot(length(L) == 1)
    stopifnot(is.numeric(L) | is.integer(L))
    stopifnot(is.character(by))
    stopifnot(by %in% c('auc', 'mapped_reads'))
    stopifnot(is.logical(factor_only))
    stopifnot(length(factor_only) == 1)
    stopifnot(is.logical(round))
    stopifnot(length(round) == 1)
    
    ## Check RSE details
    counts <- assay(rse, 1)
    stopifnot(ncol(counts) == nrow(colData(rse)))
    if(by == 'auc'){
        stopifnot('auc' %in% colnames(colData(rse)))
    } else if (by == 'mapped_reads') {
        stopifnot(all(c('avg_read_length', 'mapped_read_count',
            'paired_end') %in% colnames(colData(rse))))
    }
    
    ## Scale counts
    if(by == 'auc') {
        # L cancels out:
        # have to multiply by L to get the desired library size,
        # but then divide by L to take into account the read length since the
        # raw counts are the sum of base-level coverage.
        scaleFactor <- targetSize / colData(rse)$auc
    } else if (by == 'mapped_reads') {
        scaleFactor <- targetSize * L * ifelse(colData(rse)$paired_end, 2, 1)^2 / (colData(rse)$mapped_read_count * colData(rse)$avg_read_length^2)
    }
    
    if(factor_only) {
        names(scaleFactor) <- rownames(colData(rse))
        return(scaleFactor)
    } else {
        scaleMat <- matrix(rep(scaleFactor, each = nrow(counts)),
            ncol = ncol(counts))
        scaledCounts <- counts * scaleMat
        if(round) scaledCounts <- round(scaledCounts, 0)
        assay(rse, 1) <- scaledCounts
        return(rse)
    }
}
