#' Compute an RPKM matrix based on a RangedSummarizedExperiment object
#'
#' For some analyses you might be interested in transforming the counts into
#' RPKMs which you can do with this function.
#'
#' @param rse A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object as downloaded with [download_study].
#' @param length_var A length 1 character vector with the column name from
#' `rowData(rse)` that has the coding length. For gene level objects
#' from recount this is `bp_length`. If `NULL`, then it will use
#' `width(rowRanges(rse))` which should be used for exon RSEs.
#' @param mapped_var A length 1 character vector with the column name from
#' `colData(rse)` that has the number of reads mapped. For recount RSE
#' object this would be `mapped_read_count`. If `NULL` (default)
#' then it will use the column sums of the counts matrix. The results are
#' different because not all mapped reads are mapped to exonic segments of the
#' genome.
#'
#' @return A matrix with the RPKM values.
#'
#' @details For gene RSE objects, you will want to specify the `length_var`
#' because otherwise you will be adjusting for the total gene length instead
#' of the total exonic sequence length of the gene.
#'
#' @seealso [scale_counts]
#'
#' @author Andrew Jaffe, Leonardo Collado-Torres
#' @export
#' @import SummarizedExperiment
#'
#' @examples
#'
#' ## get RPKM matrix
#' rpkm <- getRPKM(rse_gene_SRP009615)
#'
#' ## You can also get an RPKM matrix after running scale_counts()
#' ## with similar RPKM values
#' rpkm2 <- getRPKM(scale_counts(rse_gene_SRP009615))
#' rpkm3 <- getRPKM(scale_counts(rse_gene_SRP009615, by = 'mapped_reads'))
#'
#' summary(rpkm - rpkm2)
#' summary(rpkm - rpkm3)
#' summary(rpkm2 - rpkm3)
#'

getRPKM <- function(rse, length_var = 'bp_length', mapped_var = NULL) {
    mapped <- if(!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)

    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)

    len <- if(!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))

    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    assays(rse)$counts / (wid/1000) / (bg/1e6)
}
