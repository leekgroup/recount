context('Annotation')

library('GenomicRanges')

genes <- reproduce_ranges()
exons <- reproduce_ranges('exon')

test_that('Annotation length', {
    expect_equal(length(genes), length(exons))
    expect_gte(length(genes), length(recount_genes))
    expect_gte(length(exons), length(recount_exons))
})
