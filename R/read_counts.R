#' Compute read counts
#'
#' As described in the recount workflow, the counts provided by the recount2
#' project are base-pair counts. You can scale them using \link{scale_counts}
#' or compute the read counts using the area under coverage information (AUC).
#' We use the AUC because Rail-RNA soft clips some reads.
#'
#' @param rse A \link[SummarizedExperiment]{RangedSummarizedExperiment-class} 
#' object as downloaded with \link{download_study}.
#' @param use_paired_end A logical vector. When \code{TRUE} it uses the
#' paired-end flag (\code{colData(rse)$paired_end}) to determine whether
#' the sample is paired-end or not. If it's paired-end, then it divides the
#' counts by 2 to return paired-end read counts instead of fragment counts.
#' When \code{FALSE}, this information is ignored.
#' @param round Whether to round the counts to integers or not.
#'
#' @return Returns a
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object with
#' the read counts.
#'
#' @seealso \link{scale_counts}
#'
#' @details
#' Check the recount workflow for details about the counts provided by
#' the recount2 project.
#' Note that this function should not be used after \link{scale_counts} or it 
#' will return non-sensical results.
#'
#' @references
#' Collado-Torres L, Nellore A and Jaffe AE. recount workflow: Accessing over 
#' 70,000 human RNA-seq samples with Bioconductor [version 1; referees: 1 
#' approved, 2 approved with reservations]. F1000Research 2017, 6:1558 
#' (doi: 10.12688/f1000research.12223.1)
#'
#' @author Leonardo Collado-Torres
#' @export
#' @import SummarizedExperiment
#'
#' @examples
#'
#' ## Difference between read counts and reads downloaded by Rail-RNA
#' colSums(assays(read_counts(rse_gene_SRP009615))$counts) / 1e6 - 
#'     colData(rse_gene_SRP009615)$reads_downloaded /1e6
#'
#' ## Paired-end read counts vs fragment counts (single-end counts)
#' download_study('DRP000499')
#' load('DRP000499/rse_gene.Rdata')
#' colSums(assays(read_counts(rse_gene, use_paired_end = FALSE))$counts) /
#'     colSums(assays(read_counts(rse_gene))$counts)
#'
#' ## Difference between paired-end read counts vs paired-end reads downloaded
#' colSums(assays(read_counts(rse_gene))$counts) / 1e6 - 
#'     colData(rse_gene)$reads_downloaded / 1e6 /2
#'

read_counts <- function(rse, use_paired_end = TRUE, round = FALSE) {
    
    if(use_paired_end) {
        counts <- t(t(assays(rse)$counts) /
            (colData(rse)$auc / (colData(rse)$mapped_read_count /
                ifelse(colData(rse)$paired_end, 2, 1))))
    } else {
        counts <- t(t(assays(rse)$counts) / (colData(rse)$auc /
            colData(rse)$mapped_read_count ))
    }
    if(round) counts <- round(counts, 0)
    assays(rse)$counts <- counts
    return(rse)
}
