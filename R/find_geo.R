#' Find the GEO accession id for a given SRA run
#'
#' Given a SRA run id, this function will retrieve the GEO accession id
#' (starting with GSM) if it's available. Otherwise it will return \code{NA}.
#'
#' @param run A character vector of length 1 with the SRA run accession id.
#' @param verbose Whether to print a message for the run. Useful when looping
#' over a larger number of SRA run ids.
#' @param sleep The number of seconds (or fraction) to wait before downloading
#' data using \link[GEOquery]{getGEO}. This is important if you are looking over
#' \code{geo_info()} given the constraints published at
#' http://www.ncbi.nlm.nih.gov/books/NBK25497/.
#'
#' @return The GEO accession id for the corresponding sample.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @details Although the phenotype information already includes the GEO
#' accession ids, not all projects had GEO entries at the time these tables
#' were created. This function will then be useful to check if there is a GEO
#' accession id for a given sample (run). If there is, you can then retrieve
#' the information using \link{geo_info}.
#'
#' @examples
#' ## Find the GEO accession id for for SRX110461
#' find_geo('SRX110461')
#'

find_geo <- function(run, verbose = FALSE, sleep = 1/2) {
    ## Check inputs
    stopifnot(is.character(run) & length(run) == 1)
    
    if(verbose) message(paste(Sys.time(), 'finding GEO accession id for SRA run', run))
    Sys.sleep(sleep)
    
    .load_install('XML')
    html <- XML::htmlTreeParse(paste0(
    'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=',
        run), useInternalNodes = TRUE)
    if(is.null(html)) return(NA)
    id <- XML::xpathSApply(html, '/html/body/esearchresult/idlist/id',
        XML::xmlValue)
    
    if(length(id) == 0) return(NA)
    Sys.sleep(sleep)
    
    html2 <- XML::htmlTreeParse(paste0(
        'http://www.ncbi.nlm.nih.gov/gds?LinkName=sra_gds&from_uid=', id),
        useInternalNodes = TRUE)
    if(is.null(html2)) return(NA)
    
    res <- XML::xpathSApply(html2, "//div[@class='resc']//dd", XML::xmlValue)
    gsm <- res[grep('GSM', res)]
    if(length(gsm) == 0) return(NA)
    return(gsm)
}
