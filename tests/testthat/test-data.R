context('Download data and scale')

library('SummarizedExperiment')

## Download the example data included in the package
url <- download_study('SRP009615', outdir = tempdir())

## Load the data
load(file.path(tempdir(), 'rse_gene.Rdata'))

## Compare the data
test_that('Example RSE', {
    expect_equivalent(rse_gene, rse_gene_SRP009615)
})

## Temporary test for downloading
test_that('Download URLs', {
    expect_equal(download_study('DRP000366', type = 'mean', download = FALSE),
        'http://duffel.rail.bio/recount/DRP000366/bw/mean_DRP000366.bw')
    expect_equal(download_study('DRP000366', type = 'samples',
        download = FALSE),
        'http://duffel.rail.bio/recount/DRP000366/bw/DRR000897.bw')
})

scaleFac <- scale_counts(rse_gene_SRP009615, factor_only = TRUE)
scaleFac_mapped <- scale_counts(rse_gene_SRP009615, by = 'mapped_reads', 
    factor_only = TRUE)
rse <- scale_counts(rse_gene_SRP009615, round = FALSE)

test_that('Scaling', {
    expect_equal(round(head(scaleFac), 2), c('SRR387777' = 0.04,
        'SRR387778' = 0.03, 'SRR387779' = 0.03, 'SRR387780' = 0.04,
        'SRR389077' = 0.04, 'SRR389078' = 0.04))
    expect_gt(scaleFac_mapped[1], scaleFac[1])
    expect_gt(scaleFac_mapped[2], scaleFac[2])
    expect_gt(scaleFac_mapped[3], scaleFac[3])
    expect_equal(assay(rse, 1) / matrix(rep(scaleFac, each = 23779), ncol = 12), assay(rse_gene_SRP009615, 1))
})
