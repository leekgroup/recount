context('Data')

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

test_that('Scaling', {
    expect_equal(round(head(scale_counts(rse_gene_SRP009615,
        factor_only = TRUE)), 2), c('SRR387777' = 0.04, 'SRR387778' = 0.03,
        'SRR387779' = 0.03, 'SRR387780' = 0.04, 'SRR389077' = 0.04,
        'SRR389078' = 0.04))
})
