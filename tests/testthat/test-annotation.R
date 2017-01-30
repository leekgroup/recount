context('Annotation')

library('GenomicRanges')

ranges_both <- reproduce_ranges('both', db = 'Gencode.v25')
exons <- ranges_both$exon
genes <- ranges_both$gene
exons_ensembl <- reproduce_ranges('exon', db = 'EnsDb.Hsapiens.v79')

test_that('Annotation length', {
    expect_equal(length(exons), length(genes))
    expect_gte(length(genes), length(recount_genes))
    expect_gte(length(exons), length(recount_exons))
    expect_gte(length(exons_ensembl), 65774)
    expect_equal('ENSG00000000003' %in% names(exons_ensembl), TRUE)
    expect_equal(length(genome(recount_exons)), 25)
    expect_equal('GRCh38' %in% genome(exons_ensembl), TRUE)
    expect_equal(colnames(mcols(recount_genes)), c('gene_id', 'bp_length', 'symbol'))
    expect_equal(colnames(mcols(genes)), c('gene_id', 'bp_length', 'symbol'))
})
