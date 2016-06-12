#' Extract information from GEO for a given sample
#'
#' This function uses GEOquery to extract information for a given sample. The
#' GEO accession ids for the sample can be found in the study phenotype table.
#'
#' @return Returns a \link[S4Vectors]{DataFrame-class} with the information 
#' from GEO available for the given sample.
#'
#' @param geoid A character vector of length 1 with the GEO accession id for
#' a given sample.
#' @param verbose If \code{TRUE}, messages from \link[GEOquery]{getGEO} and the
#' \code{geoid} will be shown. Otherwise they are suppressed.
#' @param sleep The number of seconds (or fraction) to wait before downloading
#' data using \link[GEOquery]{getGEO}. This is important if you are looking over
#' \code{geo_info()} given the constraints published at
#' http://www.ncbi.nlm.nih.gov/books/NBK25497/.
#' @param getGPL This argument is passed to \link[GEOquery]{getGEO} and is set
#' to \code{FALSE} by default to speed up the process.
#' @param ... Additional arguments passed to \link[GEOquery]{getGEO}. For
#' example, you might want to specify the \code{destdir} argument.
#'
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#'
#' @examples
#' geo_info('GSM836270')
#'

geo_info <- function(geoid, verbose = FALSE, sleep = 1/2, getGPL = FALSE, ...) {
    if(is.na(geoid)) return(NULL)
    
    ## Check inputs
    stopifnot(is.character(geoid) & length(geoid) == 1)
    
    ## Load GEOquery
    .load_install('GEOquery')
    .load_install('IRanges')
    
    if(verbose) message(paste(Sys.time(),
        'finding GEO information for GEO accession id', geoid))
    
    Sys.sleep(sleep)
    
    ## Get data from GEO
    if(verbose) {
        geo <- GEOquery::getGEO(geoid, getGPL = getGPL, ...)
    } else {
        geo <- suppressMessages(GEOquery::getGEO(geoid, getGPL = getGPL, ...))
    }
	
    
    ## Extract the header information
	result <- geo@header
    
    ## Function for cleaning
    clean_geo <- function(pattern, varname, res) {
    	charIndex <- grep(pattern, names(res))
    	if(length(charIndex) > 0) {
    		res <- c(res,
                IRanges::CharacterList(unlist(unname(result[charIndex]))))
            names(res)[length(res)] <- varname
    		res <- res[-charIndex]
    	}
        return(res)
    }
    
    ## Clean up the header information
    df <- data.frame(
        pattern = c('characteristics_ch1', 'data_processing', 'contact_',
            'extract_', 'library_', 'relation', 'series_',
            'supplementary_file_'),
        varname = c('characteristics', 'data_processing', 'contact', 'extract',
            'library', 'relation', 'series', 'supplementary_file'),
        stringsAsFactors = FALSE
    )        
    for(i in seq_len(nrow(df))) result <- clean_geo(df$pattern[i],
        df$varname[i], result)
    
    ## Make sure they are all length 1
    if(any(S4Vectors::elementNROWS(result) > 1)) {
        for(i in which(S4Vectors::elementNROWS(result) > 1)) result[i] <- IRanges::CharacterList(unlist(unname(result[i])))
    }
    
    ## Finish
	return(S4Vectors::DataFrame(result))
}
