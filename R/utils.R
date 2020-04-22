#' Attempt to load the namespace of one or more packages
#'
#' This function uses requireNamespace to try to load one or more packages.
#' If a package is missing, it will suggest how to install it via Bioconductor
#' before quitting.
#'
#' @param pkg A character vector with the names of the packages to check.
#'
#' @details Taken from the regionReport package
#' Updated after feedback from Marcel Ramos at
#' <https://github.com/leekgroup/recount/issues/14>
#'
#' @author Leonardo Collado-Torres
#' @noRd
#'

.load_check <- function(pkg) {
    ## Check the input
    stopifnot(is.character(pkg))
    ## Try to load the packages
    works <- sapply(pkg, requireNamespace, quietly = TRUE)

    ## If some failed, then give a useful error before quitting
    if (any(!works)) {
        x <- pkg[!works]
        stop(paste0(
            Sys.time(),
            " Package",
            ifelse(length(x) > 1, "s", ""),
            " ",
            paste(x, collapse = ", "),
            " ",
            ifelse(length(x) > 1, "are", "is"),
            " missing. Install ",
            ifelse(length(x) > 1, "them", "it"),
            " with BiocManager::install(",
            ifelse(length(x) > 1, "c(", ""),
            '"',
            paste(x, collapse = '", "'),
            '")',
            ifelse(length(x) > 1, ")", "")
        ))
    }
    return(invisible(NULL))
}
