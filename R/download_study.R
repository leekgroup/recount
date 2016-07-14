#' Download data for a given SRA study id from the recount project
#'
#' Download the gene or exon level 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment-class} objects 
#' provided by the recount project. Alternatively download the counts, metadata 
#' or file information for a given SRA study id. You can also download the 
#' sample bigWig files or the mean coverage bigWig file.
#'
#' @param project A character vector with one SRA study id.
#' @param type Specifies which files to download. The options are:
#' \describe{
#'     \item{rse-gene}{ the gene-level 
#'     \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object in 
#'     a file named rse_gene.Rdata.}
#'     \item{rse-exon}{ the exon-level 
#'     \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object in 
#'     a file named rse_exon.Rdata.}
#'     \item{counts-gene}{ the gene-level counts in a tsv file named
#'     counts_gene.tsv.gz.}
#'     \item{counts-exon}{ the exon-level counts in a tsv file named
#'     counts_exon.tsv.gz.}
#'     \item{phenotype}{ the phenotype data for the study in a tsv file named
#'     \code{project}.tsv.}
#'     \item{files-info}{ the files information for the given study (including
#'     md5sum hashes) in a tsv file named files_info.tsv.}
#'     \item{samples}{ one bigWig file per sample in the study.}
#'     \item{mean}{ one mean bigWig file for the samples in the study, 
#'     with each sample normalized to a 40 million 100 bp library using the
#'     total coverage sum (area under the coverage curve, AUC) for the given 
#'     sample.}
#'     \item{all}{ Downloads all the above types. Note that it might take some
#'     time if the project has many samples. When using \code{type = 'all'} a
#'     small delay will be added before each download request to avoid
#'     request issues.}
#' }
#' @param outdir The destination directory for the downloaded file(s).
#' @param download Whether to download the files or just get the download urls.
#' @param ... Additional arguments passed to \link[downloader]{download}.
#'
#' @return Returns invisibly the URL(s) for the files that were downloaded.
#'
#' @details Check \url{http://stackoverflow.com/a/34383991} if you need to find
#' the effective URLs. For example,
#' \url{http://duffel.rail.bio/recount/DRP000366/bw/mean_DRP000366.bw} points to
#' a temporary link from Amazon Cloud Drive.
#'
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom stats runif
#' @import downloader
#'
#' @examples
#' ## Find the URL to download the RangedSummarizedExperiment for the
#' ## Geuvadis consortium study.
#' url <- download_study('ERP001942', download = FALSE)
#' 
#' ## See the actual URL
#' url
#' 
#' ## Download the example data included in the package for study SRP009615
#'
#' url2 <- download_study('SRP009615')
#' url2
#'
#' ## Load the data
#' load(file.path('SRP009615', 'rse_gene.Rdata'))
#'
#' ## Compare the data
#' library('testthat')
#' expect_equivalent(rse_gene, rse_gene_SRP009615)
#'

download_study <- function(project, type = 'rse-gene', outdir = project,
    download = TRUE, ...) {
    ## Check inputs
    stopifnot(is.character(project) & length(project) == 1)
    stopifnot(type %in% c('rse-gene', 'rse-exon', 'counts-gene', 'counts-exon', 
        'phenotype', 'files-info', 'samples', 'mean', 'all'))
    stopifnot(length(type) == 1)
    stopifnot(is.logical(download) & length(download) == 1)
    
    ## Use table from the package
    url_table <- recount::recount_url
    
    ## Subset url data
    url_table <- url_table[url_table$project == project, ]
    stopifnot(nrow(url_table) > 0)
    
    ## If all, download each type individually
    if(type == 'all') {
        urls <- sapply(c('rse-gene', 'rse-exon', 'counts-gene', 'counts-exon', 
        'phenotype', 'files-info', 'samples', 'mean'), function(file_type) {
            Sys.sleep(round(runif(1, 2, 5), 0))
            download_study(project = project, type = file_type,
                outdir = outdir, download = download, ...)
        })
        return(invisible(urls))
    }
    
    ## Create output directory if needed
    if(download) {
        dir.create(outdir, showWarnings = FALSE)
        if(type %in% c('samples', 'mean')) {
            ## Save bigwigs on their own folder
            outdir <- file.path(outdir, 'bw')
            dir.create(outdir, showWarnings = FALSE)
        }
    }

    ## Select files to download and download them
    if(type != 'samples') {
        filename <- switch(type,
            'rse-gene' = 'rse_gene.Rdata',
            'rse-exon' = 'rse_exon.Rdata',
            'counts-gene' = 'counts_gene.tsv.gz',
            'counts-exon' = 'counts_exon.tsv.gz',
            phenotype = paste0(project, '.tsv'),
            'files-info' = 'files_info.tsv',
            mean = paste0('mean_', project, '.bw')
        )
        url <- url_table$url[url_table$file_name == filename]
        if(download) {
            message(paste(Sys.time(), 'downloading file', filename, 'to', 
                outdir))
            xx <- downloader::download(url, destfile = file.path(outdir,
                filename), mode = 'wb', ...)
        }
    } else if(type == 'samples') {
        url_table <- url_table[url_table$file_name != paste0('mean_', project,
            '.bw'), ]
        sample_urls <- url_table[grep('[.]bw$', url_table$file_name), ]
        url <- sample_urls$url
        if(download) {
            xx <- sapply(seq_len(nrow(sample_urls)), function(i, ...) {
                message(paste(Sys.time(), 'downloading file',
                    sample_urls$file_name[i], 'to', outdir))
                downloader::download(sample_urls$url[i], 
                    destfile = file.path(outdir, sample_urls$file_name[i]),
                    mode = 'wb', ...)
            }, ...)
        }
    }
    
    ## Return the actual url(s)
    return(invisible(url))
}
