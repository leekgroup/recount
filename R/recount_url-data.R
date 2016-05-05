#' Files and URLs hosted by the recount project
#'
#' Files and URLs as provided by the recount project. This information is used
#' internally in \link{download_study}.
#'
#' @name recount_url
#' @docType data
#' @format  A data.frame with 4 columns.
#' \describe{
#'     \item{path }{ the original path to the file before being uploaded,}
#'     \item{file_name }{ the file name,}
#'     \item{project }{ the SRA project id,}
#'     \item{url }{ the public URL for the given file.}
#' }
#'
#' @keywords datasets
#' @seealso \link{download_study}
#' @references \url{https://lcolladotor.shinyapps.io/recount/}
NULL 
