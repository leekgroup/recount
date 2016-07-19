#' Gene annotation used in recount
#'
#' Gene annotation extracted from TxDb.Hsapiens.UCSC.hg38.knownGene used in
#' recount. It includes the sum of the width of the reduced exons which can be
#' used for normalizing the counts provided in the 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} objects.
#'
#' @name recount_genes
#' @docType data
#' @format A \link[GenomicRanges]{GRanges-class} with one range per gene. The 
#' names match their UCSC ids. The \link[GenomicRanges]{GRanges-class} has two 
#' metadata columns.
#' \describe{
#'     \item{gene_id }{  the UCSC ids, identical to the names.}
#'     \item{bp_length }{ the sum of the width of the reduced exons for that 
#'     given gene.}
#' }
#' @keywords datasets
#' @seealso \link{reproduce_ranges}, \link{recount_exons}
#' @references \url{https://jhubiostatistics.shinyapps.io/recount/}
NULL
