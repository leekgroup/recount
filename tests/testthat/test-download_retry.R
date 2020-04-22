test_that("Downloading with retries", {
    expect_error(download_retry("http://duffel.rail.bio/recount/v2/DRP000366/files_info.tsv_fake", N.TRIES = 1))
    expect_error(download_retry("http://duffel.rail.bio/recount/v2/DRP000366/files_info.tsv", N.TRIES = 0))
    expect_error(download_retry("http://duffel.rail.bio/recount/v2/DRP000366/files_info.tsv", N.TRIES = -1L))
})
