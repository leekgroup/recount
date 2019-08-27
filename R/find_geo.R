#' Find the GEO accession id for a given SRA run
#'
#' Given a SRA run id, this function will retrieve the GEO accession id
#' (starting with GSM) if it's available. Otherwise it will return `NA`.
#'
#' @param run A character vector of length 1 with the SRA run accession id.
#' @param verbose Whether to print a message for the run. Useful when looping
#' over a larger number of SRA run ids.
#' @param sleep The number of seconds (or fraction) to wait before downloading
#' data using [getGEO][GEOquery::getGEO]. This is important if you are looking over
#' `geo_info()` given the constraints published at
#' <https://www.ncbi.nlm.nih.gov/books/NBK25497/>.
#'
#' @return The GEO accession id for the corresponding sample.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import rentrez
#'
#' @details Although the phenotype information already includes the GEO
#' accession ids, not all projects had GEO entries at the time these tables
#' were created. This function will then be useful to check if there is a GEO
#' accession id for a given sample (run). If there is, you can then retrieve
#' the information using [geo_info].
#'
#' @examples
#' ## Find the GEO accession id for for SRX110461
#' find_geo('SRX110461')
#'

find_geo <- function(run, verbose = FALSE, sleep = 1/2) {
    ## Check inputs
    stopifnot(is.character(run) & length(run) == 1)
    if(run == '') return(NA)
    
    if(verbose) message(paste(Sys.time(), 'finding GEO accession id for SRA run', run))
    Sys.sleep(sleep)
    
    ## Find uid first
    uid <- rentrez::entrez_search('sra', term = run)
    if(length(uid$ids) == 0) return(NA)
    
    ## Find linking ids
    linking <- rentrez::entrez_link('sra', id = uid$ids, db = 'gds')
    if(length(linking$links$sra_gds) == 0) return(NA)
        
    ## Find GSM
    foundGSM <- FALSE
    for(i in linking$links$sra_gds) {
        gsm <- rentrez::entrez_summary(db = 'gds', i)$accession
        if(grepl('GSM', gsm)) {
            foundGSM <- TRUE
            break
        }
    }
    
    ## Finish
    if(foundGSM) {
        return(gsm)
    } else {
        return(NA)
    }
}
