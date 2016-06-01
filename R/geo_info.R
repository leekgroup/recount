#' Extract information from GEO for a given sample
#'
#' This function uses GEOquery to extract information for a given sample. The
#' GEO accession ids for the sample can be found in the study phenotype table.
#'
#' @return Returns a list with the information from GEO available for the given
#' sample.
#'
#' @param geoid A character vector of length 1 with the GEO accession id for
#' a given sample.
#' @param verbose If \code{TRUE}, messages from \link[GEOquery]{getGEO} will
#' be shown. Otherwise they are suppressed.
#'
#' @author Leonardo Collado-Torres, Andrew Jaffe
#' @export
#'
#' @examples
#' geo_info('GSM836270')
#'

geo_info <- function(geoid, verbose = FALSE) {
    
    ## For R CMD check
    CharacterList <- getGEO <- NULL
    
    ## Check inputs
    stopifnot(is.character(geoid) & length(geoid) == 1)
    
    ## Load GEOquery
    .load_install('GEOquery')
    .load_install('IRanges')
    
    ## Get data from GEO
    if(verbose) {
        geo <- getGEO(geoid, getGPL = FALSE)
    } else {
        geo <- suppressMessages(getGEO(geoid, getGPL = FALSE))
    }
	
    
    ## Extract the header information
	result <- geo@header    
    
    ## Function for cleaning
    clean_geo <- function(pattern, varname, res) {
    	charIndex <- grep(pattern, names(res))
    	if(length(charIndex) > 0) {
    		res <- c(res, CharacterList(unname(result[charIndex])))
            names(res)[length(res)] <- varname
    		res <- res[-charIndex]
    	}
        return(res)
    }
    
    ## Clean up the header information
    df <- data.frame(
        pattern = c('characteristics_ch1', 'data_processing', 'contact_', 'extract_', 'library_', 'relation_', 'series_'),
        varname = c('characteristics', 'data_processing', 'contact', 'extract', 'library', 'relation', 'series'), stringsAsFactors = FALSE
    )        
    for(i in seq_len(nrow(df))) result <- clean_geo(df$pattern[i],
        df$varname[i], result)
    
    ## Finish
	return(result)
}
