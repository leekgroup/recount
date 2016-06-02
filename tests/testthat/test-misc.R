context('Misc')

test_that('Loading packages', {
    expect_equal(recount:::.load_install('recount'), NULL)
    expect_error(recount:::.load_install('somecrazyname'))
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
})
