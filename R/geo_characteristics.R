#' Build a data.frame from GEO's charactersitics for a given sample
#'
#' This function builds a data.frame from the GEO characteristics extracted
#' for a given sample. The names of the of columns correspond to the field
#' names. For a given SRA project, this information can be combined for all
#' samples as shown in the examples section.
#'
#' @param pheno A \link[S4Vectors]{DataFrame-class} as created by
#' \link{geo_info}.
#' 
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' ## Load required library
#' library('SummarizedExperiment')
#' 
#' ## Get the GEO accession ids
#' geoids <- sapply(colData(rse_gene_SRP009615)$run[1:2], find_geo)
#' 
#' ## Get the data from GEO
#' geodata <- do.call(rbind, sapply(geoids, geo_info))
#' 
#' ## Add characteristics in a way that we can access easily later on
#' geodata <- cbind(geodata, geo_characteristics(geodata))
#' 
#' ## Explore the original characteristics and the result from 
#' ## geo_characteristics()
#' geodata[, c('characteristics', 'cells', 'shrna.expression', 'treatment')]
#'

geo_characteristics <- function(pheno) {
    ## Check inputs
    stopifnot('characteristics' %in% colnames(pheno))
    
    .load_install('S4Vectors')
    
    ## Separate information
    result <- lapply(pheno$characteristics, function(sampleinfo) {
        info <- strsplit(sampleinfo, ': ')
    
        ## Get variable names
        varNames <- sapply(info, '[[', 1)
        varNames <- make.names(tolower(varNames))
    
        ## Construct result
        res <- matrix(sapply(info, '[[', 2), nrow = 1)
        colnames(res) <- varNames
        res <- data.frame(res, stringsAsFactors = FALSE)
    
        ## Finish
        return(res)
    })
    
    ## Finish
    result <- do.call(rbind, result) 
    return(result)
}