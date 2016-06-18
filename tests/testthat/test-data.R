context('Download data, scale, ER-level')

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
tmpdir <- file.path(tempdir(), 'SRP002001')
urls <- download_study('SRP002001', type = 'all', outdir = tmpdir)
expected_urls <- paste0('http://duffel.rail.bio/recount/SRP002001/',
    c('rse_gene.Rdata', 'rse_exon.Rdata', 'counts_gene.tsv.gz',
        'counts_exon.tsv.gz', 'SRP002001.tsv', 'files_info.tsv',
        'bw/SRR036661.bw', 'bw/mean_SRP002001.bw'))
names(expected_urls) <- c('rse-gene', 'rse-exon', 'counts-gene', 'counts-exon',
    'phenotype', 'files-info', 'samples', 'mean')

## Compute md5sum locally
localfiles <- c(list.files(tmpdir, '[.]', full.names = TRUE),
    dir(file.path(tmpdir, 'bw'), full.names = TRUE))
names(localfiles) <- c(list.files(tmpdir, '[.]'), dir(file.path(tmpdir, 'bw')))

library('tools')
md5 <- sapply(localfiles, md5sum)
names(md5) <- names(localfiles)
md5 <- md5[-which(names(md5) == 'files_info.tsv')]

## Get original md5sum
fileinfo <- read.table(file.path(tmpdir, 'files_info.tsv'), header = TRUE,
    stringsAsFactors = FALSE)

test_that('Project SRP002001', {
    expect_equal(urls, expected_urls)
    expect_equivalent(fileinfo$md5sum[match(names(md5), fileinfo$file)], md5)
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

if(.Platform$OS.type != 'windows') {
    regions <- expressed_regions('SRP002001', 'chrY', cutoff = 5)
    ## Artificially remove the mean coverage file so that the file will have to
    ## get downloaded on the first test, then it'll be present for the second
    ## test
    unlink(localfiles['mean_SRP002001.bw'])

    test_that('Expressed regions', {
        expect_equal(regions,
            expressed_regions('SRP002001', 'chrY', cutoff = 5, outdir = tmpdir))
        expect_equal(regions,
            expressed_regions('SRP002001', 'chrY', cutoff = 5, outdir = tmpdir))
    })


    coverageMatrix <- coverage_matrix('SRP002001', 'chrY', regions)
    ## Same for the phenotype data and the sample bigwig file
    unlink(localfiles['SRP002001.tsv'])
    unlink(localfiles['SRR036661.bw'])

    test_that('Coverage matrix', {
        expect_equal(coverageMatrix,
            coverage_matrix('SRP002001', 'chrY', regions, outdir = tmpdir))
        expect_equal(coverageMatrix,
            coverage_matrix('SRP002001', 'chrY', regions, outdir = tmpdir,
            chunksize = 500))
    })
}

metadata <- all_metadata()
test_that('All metadata', {
    expect_equal(nrow(metadata), 50099)
})

phenoFile <- download_study(project = 'SRP012289', type = 'phenotype',
    download = FALSE)
pheno <- read.table(phenoFile, header = TRUE, stringsAsFactors = FALSE,
    sep = '\t')
test_that('Correct phenotype information', {
    expect_equal(pheno$auc, 159080954)
})
