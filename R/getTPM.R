#' Compute a TPM matrix based on a RangedSummarizedExperiment object
#'
#' For some analyses you might be interested in transforming the counts into
#' TPMs which you can do with this function. This function uses the gene-level
#' RPKMs to derive TPM values (see Details).
#'
#'
#' @inheritParams getRPKM
#'
#' @return A matrix with the TPM values.
#'
#' @details For gene RSE objects, you will want to specify the `length_var`
#' because otherwise you will be adjusting for the total gene length instead
#' of the total exonic sequence length of the gene.
#'
#' As noted in https://support.bioconductor.org/p/124265/,
#' Sonali Arora et al computed TPMs in
#' https://www.biorxiv.org/content/10.1101/445601v2 using the formula:
#' TPM = FPKM / (sum of FPKM over all genes/transcripts) * 10^6
#'
#' Arora et al mention in their code that the formula comes from
#' https://doi.org/10.1093/bioinformatics/btp692; specifically
#' `1.1.1 Comparison to RPKM estimation` where they mention an important
#' assumption: Under the assumption of uniformly distributed reads, we note
#' that RPKM measures are estimates of ...
#'
#' There's also a blog post by Harold Pimentel explaining the relationship
#' between FPKM and TPM:
#' https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/.
#'
#' @seealso [getRPKM]
#'
#' @author Sonali Arora, Leonardo Collado-Torres
#' @export
#' @references https://www.biorxiv.org/content/10.1101/445601v2
#' https://arxiv.org/abs/1104.3889
#'
#' @examples
#'
#' ## Compute the TPM matrix from the raw gene-level base-pair counts.
#' tpm <- getTPM(rse_gene_SRP009615)
getTPM <- function(rse, length_var = "bp_length", mapped_var = NULL) {
    # From https://support.bioconductor.org/p/124265/
    rpkm <- getRPKM(rse, length_var = length_var, mapped_var = mapped_var)
    apply(rpkm, 2, function(x) {
        (x / sum(x)) * 10^6
    })
}
