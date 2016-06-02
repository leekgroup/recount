context('Misc')

test_that('Loading packages', {
    expect_equal(recount:::.load_install('recount'), NULL)
    expect_error(recount:::.load_install('somecrazyname'))
})

info <- geo_info('GSM836270')
test_that('Geo info', {
    expect_equal(info, geo_info('GSM836270', TRUE))
    expect_equal(info$geo_accession, 'GSM836270')
    expect_equal(find_geo('SRX110461'), 'GSM836270')
    expect_equal(find_geo(''), NA)
})
