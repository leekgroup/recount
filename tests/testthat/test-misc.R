context('Misc')

test_that('Loading packages', {
    expect_equal(recount:::.load_install('recount'), NULL)
    expect_error(recount:::.load_install('somecrazyname'))
})