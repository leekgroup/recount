context('Misc')

test_that('Loading packages', {
    expect_equal(recount:::.load_check('recount'), NULL)
    expect_error(recount:::.load_check('somecrazyname'))
})

info <- geo_info('GSM836270')
chars <- data.frame('cells' = 'K562', 'shrna.expression' = 'no',
    'treatment' = 'Puromycin', stringsAsFactors = FALSE)
test_that('Geo info', {
    expect_equal(info, geo_info('GSM836270', TRUE))
    expect_equal(info$geo_accession, 'GSM836270')
    expect_equal(find_geo('SRX110461'), 'GSM836270')
    expect_equal(find_geo(''), NA)
    expect_equal(geo_characteristics(info), chars)
    expect_equal(geo_info(NA), NULL)
    expect_equal(colnames(geo_characteristics(geo_info('GSM359183'))), 'characteristics')
    #expect_equal(geo_info('GSM1062236'), S4Vectors::DataFrame())
})

if(interactive()) {
    ## Open shiny app when running test in interactive mode.
    ## This will force the app to be available, which will then run the test
    ## successfully.
    ## If there is a better way to test this, please let me know!
    ## See thread at https://groups.google.com/forum/#!topic/shinyapps-users/ElMO_v1eurQ
    browseURL('https://jhubiostatistics.shinyapps.io/recount/')
    test_that('Shiny app is up', {
        expect_equal(
            RCurl::url.exists('https://jhubiostatistics.shinyapps.io/recount/'),
            TRUE
        )
    })
}

