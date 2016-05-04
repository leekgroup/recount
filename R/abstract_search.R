#' Search the abstracts from the SRA studies available via the recount project
#'
#' Given a text query, find the SRA project ids (study accession numbers) that
#' contain the text in their abstract as provided by the SRAdb Bioconductor
#' package.
#'
#' @param query A character vector with the text to search for via
#' \link[base]{grep} in the abstract info available at \link{recount_abstract}.
#' @param id_only Whether to only return the project id or to return summary
#' information for the project(s) that match the query.
#' @param ... Additional arguments passed to \link[base]{grep}.
#'
#' @return If \code{id_only = TRUE} it returns a character vector with the
#' project SRA ids (accession numbers). If \code{id_only = FALSE} it returns a
#' subset of \link{recount_abstract} for the abstracts that contained the query.
#'
#' @details Both the query and the abstracts are searched in lower case.
#'
#' For a more powerful search use the recount project website at
#' \url{https://lcolladotor.shinyapps.io/recount/}.
#'
#' @seealso \link{browse_study}, \link{recount_abstract}
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' ## Find the Geuvadis consortium project
#' project_info <- abstract_search('Geuvadis consortium')
#' 
#' ## See some summary information for this project
#' project_info

abstract_search <- function(query, id_only = FALSE, ...) {
    ## Check input
    stopifnot(is.character(query))
    stopifnot(is.logical(id_only))
    stopifnot(length(id_only) == 1)
    
    ## Use table from the package
    abstract_table <- recount::recount_abstract
    
    ## Get abstracts
    abstracts <- tolower(abstract_table$abstract)
    query <- tolower(query)
    i <- grep(query, abstracts, ...)
    
    if(id_only) {
        return(abstract_table$project[i])
    } else {
        return(abstract_table[i, ])
    }
}
