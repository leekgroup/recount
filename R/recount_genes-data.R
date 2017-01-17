#' Gene annotation used in recount
#'
#' Gene annotation extracted from Gencode v25 (GRCh38.p7) used in recount.
#' Data extraced on January 17th, 2017.
#' It includes the sum of the width of the reduced exons which can be
#' used for normalizing the counts provided in the 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} objects.
#'
#' @name recount_genes
#' @docType data
#' @format A \link[GenomicRanges]{GRanges-class} with one range per gene. The 
#' names match their Gencode v25 ids. The \link[GenomicRanges]{GRanges-class}
#' has three metadata columns.
#' \describe{
#'     \item{gene_id }{  the Gencode v25 ids, identical to the names.}
#'     \item{bp_length }{ the sum of the width of the reduced exons for that 
#'     given gene.}
#'     \item{symbol }{ a CharacterList with the corresponding gene symbols.}
#' }
#' @keywords datasets
#' @seealso \link{reproduce_ranges}, \link{recount_exons}
#' @references \url{https://jhubiostatistics.shinyapps.io/recount/}
NULL
