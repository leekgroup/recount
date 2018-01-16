context('Website metadata')

test_that('Browse studies at SRA', {
    expect_equal(length(browse_study(c('ERP001942', 'ERP009768'), FALSE)), 2)
    expect_equal(browse_study('ERP001942', FALSE), 'https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP001942')
})

test_that('Download studies', {
    expect_equal(download_study('ERP001942', download = FALSE), 'http://duffel.rail.bio/recount/v2/ERP001942/rse_gene.Rdata')
    expect_equal(download_study('ERP001942', download = FALSE, version = 1), 'http://duffel.rail.bio/recount/ERP001942/rse_gene.Rdata')
    expect_equal(download_study('ERP001942', type = 'mean', download = FALSE), 'http://duffel.rail.bio/recount/ERP001942/bw/mean_ERP001942.bw')
})


test_that('Finding abstracts', {
    expect_equal(abstract_search('Geuvadis consortium', id_only = TRUE),
        'ERP001942')
    expect_equal(abstract_search('Geuvadis consortium'),
        recount_abstract[recount_abstract$project == 'ERP001942', ])
})
