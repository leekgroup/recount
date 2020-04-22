#' Extract information from GEO for a given sample
#'
#' This function uses GEOquery to extract information for a given sample. The
#' GEO accession ids for the sample can be found in the study phenotype table.
#'
#' @return Returns a [DataFrame-class][S4Vectors::DataFrame-class] with the information
#' from GEO available for the given sample.
#'
#' @param geoid A character vector of length 1 with the GEO accession id for
#' a given sample.
#' @param verbose If `TRUE` the `geoid` will be shown.
#' @param sleep The number of seconds (or fraction) to wait before downloading
#' data using [getGEO][GEOquery::getGEO]. This is important if you are looking over
#' `geo_info()` given the constraints published at
#' <https://www.ncbi.nlm.nih.gov/books/NBK25497/>.
#' @param getGPL This argument is passed to [getGEO][GEOquery::getGEO] and is set
#' to `FALSE` by default to speed up the process.
#' @param destdir This argument is passed to [getGEO][GEOquery::getGEO].
#' @param ... Additional arguments passed to [getGEO][GEOquery::getGEO].
#'
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#'
#' @import GEOquery IRanges S4Vectors
#'
#' @examples
#' geo_info("GSM836270")
geo_info <- function(geoid, verbose = FALSE, sleep = 1 / 2, getGPL = FALSE,
    destdir = tempdir(), ...) {
    if (is.na(geoid)) {
        return(NULL)
    }

    ## Check inputs
    stopifnot(is.character(geoid) & length(geoid) == 1)

    if (verbose) {
        message(paste(
            Sys.time(),
            "finding GEO information for GEO accession id", geoid
        ))
    }

    if (!file.exists(file.path(destdir, paste0(geoid, ".soft")))) {
        Sys.sleep(sleep)
    }

    ## Get data from GEO, with 3 retries, waiting between 0 and 2 seconds in
    ## between retries
    N.TRIES <- 3L
    while (N.TRIES > 0L) {
        geo <- tryCatch(
            GEOquery::getGEO(geoid, getGPL = getGPL, destdir = destdir, ...),
            error = function(e) {
                soft <- paste0(geoid, ".soft")
                soft_file <- file.path(destdir, soft)
                if (any(grepl("private", readLines(soft_file)))) {
                    message(paste(geoid, "is currently private"))
                    return(NA)
                } else if (any(grepl("blocked", readLines(soft_file)))) {
                    warning(paste("It seems like your IP access is blocked. Please check the file", soft_file, "for more details."))
                    return(NA)
                } else {
                    return(e)
                }
            }
        )
        if (!inherits(geo, "error")) {
              break
          }
        Sys.sleep(runif(n = 1, min = 2, max = 5))
        N.TRIES <- N.TRIES - 1L
    }

    ## Return and empty DataFrame if there was an issue with getGEO()
    if (!is(geo, "GSM")) {
        return(S4Vectors::DataFrame())
    }

    ## Extract the header information
    result <- geo@header

    ## Function for cleaning
    clean_geo <- function(pattern, varname, res) {
        charIndex <- grep(pattern, names(res))
        if (length(charIndex) > 0) {
            res <- c(
                res,
                IRanges::CharacterList(unlist(unname(result[charIndex])))
            )
            names(res)[length(res)] <- varname
            res <- res[-charIndex]
        }
        return(res)
    }

    ## Clean up the header information
    df <- data.frame(
        pattern = c(
            "characteristics_ch1", "data_processing", "contact_",
            "extract_", "library_", "relation", "series_",
            "supplementary_file_"
        ),
        varname = c(
            "characteristics", "data_processing", "contact", "extract",
            "library", "relation", "series", "supplementary_file"
        ),
        stringsAsFactors = FALSE
    )
    for (i in seq_len(nrow(df))) {
        result <- clean_geo(
            df$pattern[i],
            df$varname[i], result
        )
    }

    ## Make sure they are all length 1
    if (any(S4Vectors::elementNROWS(result) > 1)) {
        for (i in which(S4Vectors::elementNROWS(result) > 1)) result[i] <- IRanges::CharacterList(unlist(unname(result[i])))
    }

    ## Finish
    return(S4Vectors::DataFrame(result))
}
