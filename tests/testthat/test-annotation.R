context('Annotation')

library('GenomicRanges')

genes <- reproduce_ranges()
exons <- reproduce_ranges('exon')

exons_ensembl <- reproduce_ranges('exon', db = 'EnsDb.Hsapiens.v79')

test_that('Annotation length', {
    expect_equal(length(genes), length(exons))
    expect_gte(length(genes), length(recount_genes))
    expect_gte(length(exons), length(recount_exons))
    expect_gte(length(exons_ensembl), 65774)
    expect_equal('ENSG00000000003' %in% names(exons_ensembl), TRUE)
    expect_equal('hg38' %in% genome(recount_exons), TRUE)
    expect_equal('hg38' %in% genome(recount_genes), TRUE)
    expect_equal('GRCh38' %in% genome(exons_ensembl), TRUE)
})
