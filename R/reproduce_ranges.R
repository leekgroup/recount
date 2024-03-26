#' Reproduce the gene or exons used in the RangedSummarizedExperiment objects
#'
#' This function reproduces the gene or exon level information used for creating
#' the [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' objects provided by recount. The annotation is based on
#' Gencode v25 with the gene-level
#' information extracted with `genes()` (see
#' [transcripts][GenomicFeatures::transcripts] with default arguments.
#'
#' @param level Either `genes` or `exon`. It specifies whether to
#' return Gene or exon level information as a
#' [GRanges-class][GenomicRanges::GRanges-class] or
#' [GRangesList-class][GenomicRanges::GRangesList-class] object respectively. The gene level
#' information contains the width of the disjoint exons for that given gene
#' which can be used to normalize the counts provided by recount.
#' Can also be `both` in which case a 2 element list with the exon and the
#' gene output is returned.
#' @param db Either `Gencode.v25` (default) or
#' `EnsDb.Hsapiens.v79`. The default option reproduces the annotation
#' used when creating recount. EnsDb.Hsapiens.v79 can be used
#' for an alternative annotation as showcased in the recount vignette.
#'
#' @return Either a [GRanges-class][GenomicRanges::GRanges-class] object like
#' [recount_genes] or a [GRangesList-class][GenomicRanges::GRangesList-class] object
#' like [recount_exons].
#'
#' @details
#'
#' For Gencode.v25, we use the comprehensive gene annotation (regions:
#' `CHR`) from <https://www.gencodegenes.org/releases/25.html>
#' (GRCh38.p7).
#'
#' Note that gene symbols have changed over time. This answer in the
#' Bioconductor support forum details how to obtain the latest gene symbol
#' mappings: <https://support.bioconductor.org/p/126148/#126173>.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import GenomicRanges
#'
#' @seealso [recount_genes], [recount_exons],
#'     <https://github.com/nellore>,
#'     <https://jhubiostatistics.shinyapps.io/recount/>
#'
#' @examples
#' \dontrun{
#' ## Reproduce gene level information
#' genes <- reproduce_ranges()
#'
#' ## Compare against recount_genes
#' length(genes)
#' length(recount_genes)
#' }
#'
reproduce_ranges <- function(level = "gene", db = "Gencode.v25") {
    ## Check input
    level <- tolower(level)
    stopifnot(level %in% c("gene", "exon", "both"))
    stopifnot(db %in% c("Gencode.v25", "EnsDb.Hsapiens.v79"))



    ## Load required packages
    .load_check(c("GenomicFeatures", "org.Hs.eg.db"))
    if (db == "Gencode.v25") {
        .load_check("txdbmaker")
        temp_gencode <- file.path(
            tempdir(),
            "gencode.v25.annotation.gff3.gz"
        )
        xx <- download_retry(
            url = paste0(
                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/",
                "Gencode_human/release_25/gencode.v25.annotation.gff3.gz"
            ),
            destfile = temp_gencode
        )
        txdb <- txdbmaker::makeTxDbFromGFF(temp_gencode,
            format = "gff3", organism = "Homo sapiens"
        )
    } else if (db == "EnsDb.Hsapiens.v79") {
        .load_check("EnsDb.Hsapiens.v79")
        txdb <- EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79
    }

    ## Get genes with default option single.strand.genes.only = TRUE
    genes <- GenomicFeatures::genes(txdb)

    ## Get Exons
    exons <- GenomicFeatures::exonsBy(txdb, by = "gene")

    ## Keep only exons for gene ids that we selected previously
    if (!all(names(exons) %in% names(genes))) {
        warning("Dropping exons with gene ids not present in the gene list")
        exons <- exons[names(exons) %in% names(genes)]
    }

    ## Disjoin exons by gene so the exons won't be overlapping each other inside a gene
    exons <- GenomicRanges::disjoin(exons)

    if (level == "exon") {
        return(exons)
    }

    ## For 'gene' or 'both', continue by:
    ## * adding length of disjoint exons by gene
    genes$bp_length <- sum(GenomicRanges::width(exons))

    ## * adding gene symbol
    if (db == "Gencode.v25") {
        gene_info <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
            gsub("\\..*", "", names(genes)), "SYMBOL", "ENSEMBL",
            multiVals = "CharacterList"
        )
    } else if (db == "EnsDb.Hsapiens.v79") {
        gene_info <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
            names(genes), "SYMBOL", "ENSEMBL",
            multiVals = "CharacterList"
        )
    }
    genes$symbol <- gene_info

    if (level == "gene") {
        return(genes)
    } else if (level == "both") {
        return(list("exon" = exons, "gene" = genes))
    }
}
