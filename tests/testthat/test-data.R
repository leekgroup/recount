context('Download data and scale')

library('SummarizedExperiment')

## Download the example data included in the package
url <- download_study('SRP009615', outdir = file.path(tempdir(), 'SRP009615'))

## Load the data
load(file.path(tempdir(), 'SRP009615', 'rse_gene.Rdata'))

## Compare the data
test_that('Example RSE', {
    expect_equivalent(rse_gene, rse_gene_SRP009615)
})

## Temporary download URLs
test_that('Download URLs', {
    expect_equal(download_study('DRP000366', type = 'mean', download = FALSE),
        'http://duffel.rail.bio/recount/DRP000366/bw/mean_DRP000366.bw')
    expect_equal(download_study('DRP000366', type = 'samples',
        download = FALSE),
        'http://duffel.rail.bio/recount/DRP000366/bw/DRR000897.bw')
})

## Test downloading a small project entirely
urls <- download_study('SRP002001', type = 'all', outdir = file.path(tempdir(),
    'SRP002001'))
expected_urls <- paste0('http://duffel.rail.bio/recount/SRP002001/',
    c('rse_gene.Rdata', 'rse_exon.Rdata', 'counts_gene.tsv.gz',
        'counts_exon.tsv.gz', 'SRP002001.tsv', 'files_info.tsv',
        'bw/SRR036661.bw', 'bw/mean_SRP002001.bw'))
names(expected_urls) <- c('rse-gene', 'rse-exon', 'counts-gene', 'counts-exon',
    'phenotype', 'files-info', 'samples', 'mean')

test_that('Project SRP002001', {
    expect_equal(urls, expected_urls)
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
    expect_equal(assay(rse, 1) / matrix(rep(scaleFac, each = 23779),
        ncol = 12), assay(rse_gene_SRP009615, 1))
})
