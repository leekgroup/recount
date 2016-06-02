#' Find the GEO accession id for a given SRA run
#'
#' Given a SRA run id, this function will retrieve the GEO accession id
#' (starting with GSM) if it's available. Otherwise it will return \code{NA}.
#'
#' @param run A character vector of length 1 with the SRA run accession id.
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

find_geo <- function(run) {
    
    ## For R CMD check
    htmlTreeParse <- xpathSApply <- NULL
    ## Can't do the same to xmlValue
    
    ## Check inputs
    stopifnot(is.character(run) & length(run) == 1)
    
    .load_install('XML')
    html <- htmlTreeParse(paste0(
    'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=',
        run), useInternalNodes = TRUE)
    id <- xpathSApply(html, '/html/body/esearchresult/idlist/id', xmlValue)
    
    if(length(id) == 0) return(NA)
    
    html2 <- htmlTreeParse(paste0(
        'http://www.ncbi.nlm.nih.gov/gds?LinkName=sra_gds&from_uid=', id),
        useInternalNodes = TRUE)
    
    res <- xpathSApply(html2, "//div[@class='resc']//dd", xmlValue)
    gsm <- res[grep('GSM', res)]
    if(length(gsm) == 0) return(NA)
    return(gsm)
}
