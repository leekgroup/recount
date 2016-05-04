#' Open a SRA study id in the SRA website
#'
#' Given a SRA study id get the url to browse the study using the SRA website.
#'
#' @param project A character vector with at least one SRA study id.
#' @param browse Whether to open the resulting URL in the browser.
#'
#' @return Returns invisibly the URL for exploring the study in the SRA website.
#'
#' @author Leonardo Collado-Torres
#' @export
#' @seealso \link{abstract_search}
#'
#' @examples
#' ## Find the Geuvadis consortium project
#' id <- abstract_search('Geuvadis consortium', id_only = TRUE)
#' id
#'
#' ## Explore the Geuvadis consortium project
#' url <- browse_study(id)
#' 
#' ## See the actual URL
#' url

browse_study <- function(project, browse = interactive()) {
    ## Check inputs
    stopifnot(is.logical(browse) & length(browse == 1))
    stopifnot(is.character(project))
    
    ## Construct url
    url <- paste0('http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=', project)
    
    ## Finish
    if(browse) sapply(url, browseURL)
    return(invisible(url))
}
