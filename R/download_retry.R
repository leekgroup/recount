#' Retry multiple times to download a file
#'
#' This function is based on the Bioconductor guidelines for querying data from
#' the web at <http://bioconductor.org/developers/how-to/web-query/>. It
#' will run [download][downloader::download] a set of `N.TRIES` times before
#' giving up. We implemented this function to reduce the number of Bioconductor
#' build errors due to the occassional errors from our data hosting server.
#'
#' @param url The URL to download. Passed to [download][downloader::download].
#' @param destfile The destination file. Defaults to the base name of the URL.
#' Passed to [download][downloader::download].
#' @param mode Mode for writing the file. The default `wb` is used for
#' binary files. This value is passed to [download][downloader::download] which
#' passes it to [download.file][utils::download.file].
#' @param N.TRIES The number of download attempts before giving up; default: 3.
#' Should be an integer of length one with a value greater than 0.
#' @param ... Additional arguments passed to [download][downloader::download].
#'
#' @return An invisible integer code as specified in
#' [download.file][utils::download.file].
#' @export
#'
#' @examples
#'
#' ## Download the first files_info.tsv file (version 1)
#' download_retry(
#'     recount_url$url[which(recount_url$file_name == 'files_info.tsv')[1]]
#' )
#'
download_retry <- function(url, destfile = basename(url), mode = 'wb',
    N.TRIES = 3L, ...) {
    ## Based on http://bioconductor.org/developers/how-to/web-query/
    ## and downloader::download()

    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
    stopifnot(N.TRIES > 0L)

    while (N.TRIES > 0L) {
        result <- tryCatch(downloader::download(
            url = url, destfile = destfile, mode = mode, ...), error=identity)
        if (!inherits(result, "error"))
            break
        ## Wait between 0 and 2 seconds between retries
        Sys.sleep(runif(n = 1, min = 2, max = 5))
        N.TRIES <- N.TRIES - 1L
    }

    if (N.TRIES == 0L) {
        stop("'download_retry()' failed:",
             "\n  URL: ", url,
             "\n  error: ", conditionMessage(result))
    }

    invisible(result)
}
