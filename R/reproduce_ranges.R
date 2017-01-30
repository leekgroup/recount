#' Reproduce the gene or exons used in the RangedSummarizedExperiment objects
#'
#' This function reproduces the gene or exon level information used for creating
#' the \link[SummarizedExperiment]{RangedSummarizedExperiment-class}
#' objects provided by recount. The annotation is based on
#' Gencode v25 with the gene-level 
#' information extracted with \code{genes()} (see 
#' \link[GenomicFeatures]{transcripts} with default arguments.
#' 
#' @param level Either \code{genes} or \code{exon}. It specifies whether to
#' return Gene or exon level information as a 
#' \link[GenomicRanges]{GRanges-class} or 
#' \link[GenomicRanges]{GRangesList-class} object respectively. The gene level
#' information contains the width of the reduced exons for that given gene
#' which can be used to normalize the counts provided by recount.
#' Can also be \code{both} in which case a 2 element list with the exon and the
#' gene output is returned.
#' @param db Either \code{Gencode.v25} (default) or
#' \code{EnsDb.Hsapiens.v79}. The default option reproduces the annotation
#' used when creating recount. EnsDb.Hsapiens.v79 can be used
#' for an alternative annotation as showcased in the recount vignette.
#'
#' @return Either a \link[GenomicRanges]{GRanges-class} object like 
#' \link{recount_genes} or a \link[GenomicRanges]{GRangesList-class} object 
#' like \link{recount_exons}.
#'
#' @details
#' 
#' For Gencode.v25, we use the comprehensive gene annotation (regions: 
#' \code{CHR}) from \url{https://www.gencodegenes.org/releases/25.html}
#' (GRCh38.p7).
#' 
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import GenomicRanges
#' @import GenomicFeatures
#'
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

reproduce_ranges <- function(level = 'gene', db = 'Gencode.v25') {
    ## Check input
    level <- tolower(level)
    stopifnot(level %in% c('gene', 'exon', 'both'))
    stopifnot(db %in% c('Gencode.v25', 'EnsDb.Hsapiens.v79'))
    
    
    
    ## Load required packages
    .load_install('GenomicFeatures')
    if (db == 'Gencode.v25') {
            txdb <- GenomicFeatures::makeTxDbFromGFF('ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gff3.gz',
            format = 'gff3', organism = 'Homo sapiens')
    } else if(db == 'EnsDb.Hsapiens.v79') {
        .load_install('EnsDb.Hsapiens.v79')
        txdb <- EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79
    }

    ## Get genes with default option single.strand.genes.only = TRUE
    genes <- GenomicFeatures::genes(txdb)

    ## Get Exons
    exons <- GenomicFeatures::exonsBy(txdb, by = 'gene')

    ## Keep only exons for gene ids that we selected previously
    if(!all(names(exons) %in% names(genes))) {
        warning('Dropping exons with gene ids not present in the gene list')
        exons <- exons[names(exons) %in% names(genes)]
    }
        
    ## Reduce exons by gene so the exons won't be overlapping each other inside a gene
    exons <- GenomicRanges::reduce(exons)

    if(level == 'exon') return(exons)
        
    ## For 'gene' or 'both', continue by:
    ## * adding length of reduced exons by gene
    genes$bp_length <- sum(GenomicRanges::width(exons))
    
    ## * adding gene symbol
    .load_install('org.Hs.eg.db')
    if(db == 'Gencode.v25') {
        gene_info <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
            gsub('\\..*', '', names(genes)), 'SYMBOL', 'ENSEMBL',
            multiVals = 'CharacterList')
    } else if(db == 'EnsDb.Hsapiens.v79') {
        gene_info <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
            names(genes), 'SYMBOL', 'ENSEMBL', multiVals = 'CharacterList')
    }
    genes$symbol <- gene_info
    
    if(level == 'gene') {
        return(genes)
    } else if (level == 'both') {
        return(list('exon' = exons, 'gene' = genes))
    }
}

