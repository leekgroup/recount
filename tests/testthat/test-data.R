context('Download data and scale')

library('SummarizedExperiment')

if(FALSE) {
    ## Test won't work until data is fully uploaded

## Download the example data included in the package
url <- download_study('SRP009615', outdir = tempdir())

## Load the data
load(file.path(tempdir(), 'SRP009615', 'rse_gene.Rdata'))

## Compare the data
test_that('Example RSE', {
    expect_equal(rse_gene, rse_gene_SRP009615)
})
}

## Temporary test for downloading
## Hm... somehow download.file('https://lcolladotor.shinyapps.io/recount/ucsc-knowngene-hg38-genes-bp-length.Rdata', destfile = 'x.Rdata') doesn't work
test_that('Downloading', {
    expect_error(download_study('DRP000366', type = 'mean', outdir = tempdir()))
    expect_error(download_study('DRP000366', type = 'samples', outdir = tempdir()))
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
