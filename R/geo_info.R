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
#' @param verbose If \code{TRUE} the \code{geoid} will be shown. 
#' @param sleep The number of seconds (or fraction) to wait before downloading
#' data using \link[GEOquery]{getGEO}. This is important if you are looking over
#' \code{geo_info()} given the constraints published at
#' http://www.ncbi.nlm.nih.gov/books/NBK25497/.
#' @param getGPL This argument is passed to \link[GEOquery]{getGEO} and is set
#' to \code{FALSE} by default to speed up the process.
#' @param destdir This argument is passed to \link[GEOquery]{getGEO}.
#' @param ... Additional arguments passed to \link[GEOquery]{getGEO}.
#'
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#'
#' @examples
#' geo_info('GSM836270')
#'

geo_info <- function(geoid, verbose = FALSE, sleep = 1/2, getGPL = FALSE,
    destdir = tempdir(), ...) {
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
    geo <- tryCatch(
        GEOquery::getGEO(geoid, getGPL = getGPL, destdir = destdir, ...), 
        error = function(e) {
            soft <- paste0(geoid, '.soft')
            soft_file <- file.path(destdir, soft)
            if(any(grepl('private', readLines(soft_file)))) {
                message(paste(geoid, 'is currently private'))
                return(NA)
            } else if (any(grepl('blocked', readLines(soft_file)))) {
                warning(paste('It seems like your IP access is blocked. Please check the file', soft_file, 'for more details.'))
                return(NA)
            } else {
                stop(e)
            }
        }
    )
	
    ## Return and empty DataFrame if there was an issue with getGEO()
    if(!is(geo, 'GSM')) return(S4Vectors::DataFrame())
    
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
