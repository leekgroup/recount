#' Files and URLs hosted by the recount project
#'
#' Files and URLs as provided by the recount project. This information is used
#' internally in [download_study].
#'
#' @name recount_url
#' @docType data
#' @format  A data.frame with 6 columns.
#' \describe{
#'     \item{path }{ the original path to the file before being uploaded,}
#'     \item{file_name }{ the file name,}
#'     \item{project }{ the SRA project id,}
#'     \item{version1 }{ A logical vector indicating whether the file was
#'     part of version 1 (reduced exons)}.
#'     \item{version 2}{ A logical vector indicating whether the file was
#'     updated in version 2 (disjoint exons)}. Further details in the recount
#'     website documentation tab.
#'     \item{url }{ the public URL for the given file.}
#' }
#'
#' @keywords datasets
#' @seealso [download_study]
#' @references <https://jhubiostatistics.shinyapps.io/recount/>
NULL 
