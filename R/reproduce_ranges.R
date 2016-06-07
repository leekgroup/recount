#' Reproduce the gene or exons used in the RangedSummarizedExperiment objects
#'
#' This function reproduces the gene or exon level information used for creating
#' the \link[SummarizedExperiment]{RangedSummarizedExperiment-class}
#' objects provided by recount. The annotation is based on
#' \code{TxDb.Hsapiens.UCSC.hg38.knownGene} with the gene-level 
#' information extracted with \code{genes()} (see 
#' \link[GenomicFeatures]{transcripts} with default arguments.
#' 
#' @param level Either \code{genes} or \code{exon}. It specifies whether to
#' return Gene or exon level information as a 
#' \link[GenomicRanges]{GRanges-class} or 
#' \link[GenomicRanges]{GRangesList-class} object respectively. The gene level
#' information contains the width of the reduced exons for that given gene
#' which can be used to normalize the counts provided by recount.
#'
#' @return Either a \link[GenomicRanges]{GRanges-class} object like 
#' \link{recount_genes} or a \link[GenomicRanges]{GRangesList-class} object 
#' like \link{recount_exons}.
#'
#' @details
#' Note that the information used in recount was originally calculated with
#' \url{https://github.com/leekgroup/recount-website/blob/master/genes/create_bed.R}
#' with the versions as noted in
#' \url{https://github.com/leekgroup/recount-website/blob/master/genes/logs/genes-hg38-ucsc.o1192789}.
#' The actual package tarball for TxDb.Hsapiens.UCSC.hg38.knownGene we used is 
#' available at \url{https://jhubiostatistics.shinyapps.io/recount/TxDb.Hsapiens.UCSC.hg38.knownGene_3.1.3.tar.gz}.
#' 
#' @author Leonardo Collado-Torres
#' @export
#' @seealso \link{recount_genes}, \link{recount_exons},
#'     \url{https://github.com/nellore},
#'     \url{https://jhubiostatistics.shinyapps.io/recount/}
#' 
#' @examples
#'
#' ## Reproduce gene level information
#' genes <- reproduce_ranges()
#' 
#' \dontrun{
#' ## Compare against recount_genes
#' length(genes)
#' length(recount_genes)
#' }
#'

reproduce_ranges <- function(level = 'gene') {
    ## Check input
    stopifnot(level %in% c('gene', 'exon'))
    
    ## Load required packages
    .load_install('GenomicFeatures')
    .load_install('GenomicRanges')
    .load_install('TxDb.Hsapiens.UCSC.hg38.knownGene')

    ## Get genes with default option single.strand.genes.only = TRUE
    genes <- GenomicFeatures::genes(
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)

    ## Get Exons
    exons <- GenomicFeatures::exonsBy(
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
        by = 'gene')

    ## Keep only exons for gene ids that we selected previously
    exons <- exons[names(exons) %in% names(genes)]
    
    ## Reduce exons by gene so the exons won't be overlapping each other inside a gene
    exons <- GenomicRanges::reduce(exons)
    
    if(level == 'exon') {
        return(exons)
    } else if(level == 'gene') {
        ## Add length of reduced exons by gene
        genes$bp_length <- sum(GenomicRanges::width(exons))
        
        return(genes)
    }
}

