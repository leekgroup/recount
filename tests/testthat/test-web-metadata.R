context('Website metadata')

test_that('Browse studies at SRA', {
    expect_equal(length(browse_study(c('ERP001942', 'ERP009768'), FALSE)), 2)
    expect_equal(browse_study('ERP001942', FALSE), 'http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP001942')
})

test_that('Download studies', {
    expect_equal(download_study('ERP001942', download = FALSE), 'https://lcolladotor.shinyapps.io/recount/ucsc-knowngene-hg38-genes-bp-length.Rdata')
})


test_that('Finding abstracts', {
    expect_equal(abstract_search('Geuvadis consortium', id_only = TRUE),
        'ERP001942')
    expect_equal(abstract_search('Geuvadis consortium'),
        recount_abstract[recount_abstract$project == 'ERP001942', ])
})
