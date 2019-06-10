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
#'     \item{rse-jx}{ the exon-exon junction level
#'     \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object in
#'     a file named rse_jx.Rdata.}
#'     \item{rse-tx}{ the transcript level
#'     \link[SummarizedExperiment]{RangedSummarizedExperiment-class} object in
#'     a file named rse_tx.RData.}
#'     \item{counts-gene}{ the gene-level counts in a tsv file named
#'     counts_gene.tsv.gz.}
#'     \item{counts-exon}{ the exon-level counts in a tsv file named
#'     counts_exon.tsv.gz.}
#'     \item{counts-jx}{ the exon-exon junction level counts in a tsv file named
#'     counts_jx.tsv.gz.}
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
#'     \item{rse-fc}{ Downloads the FANTOM-CAT/recount2 rse file described in
#'     Imada, Sanchez, et al., bioRxiv, 2019.}
#' }
#' @param outdir The destination directory for the downloaded file(s).
#' Alternatively check the \code{SciServer} section on the vignette to see
#' how to access all the recount data via a R Jupyter Notebook.
#' @param download Whether to download the files or just get the download urls.
#' @param version A single integer specifying which version of the files to
#' download. Valid options are 1 and 2, as described in
#' \url{https://jhubiostatistics.shinyapps.io/recount/} under the
#' documentation tab. Briefly, version 1 are counts based on reduced exons while
#' version 2 are based on disjoint exons. This argument mostly just matters for
#' the exon counts. Defaults to version 2 (disjoint exons).
#' Use \code{version = 1} for backward compatability with exon counts
#' prior to version 1.5.3 of the package.
#' @param ... Additional arguments passed to \link[downloader]{download}.
#'
#' @return Returns invisibly the URL(s) for the files that were downloaded.
#'
#' @details Check \url{http://stackoverflow.com/a/34383991} if you need to find
#' the effective URLs. For example,
#' \url{http://duffel.rail.bio/recount/DRP000366/bw/mean_DRP000366.bw} points to
#' a link from SciServer.
#'
#' Transcript quantifications are described in Fu et al, bioRxiv, 2018.
#' \url{https://www.biorxiv.org/content/10.1101/247346v2}
#'
#' FANTOM-CAT/recount2 quantifications are described in Imada,
#' Sanchez, et al., bioRxiv, 2019.
#' \url{https://www.biorxiv.org/content/10.1101/659490v1}
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
#' \dontrun{
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
#' }
#'

download_study <- function(project, type = 'rse-gene', outdir = project,
    download = TRUE, version = 2, ...) {
    ## Check inputs
    stopifnot(is.character(project) & length(project) == 1)
    stopifnot(version %in% c(1, 2))
    stopifnot(type %in% c(
        'rse-gene', 'rse-exon', 'rse-jx', 'rse-tx',
        'counts-gene', 'counts-exon', 'counts-jx',
        'phenotype', 'files-info', 'samples', 'mean', 'rse-fc', 'all'))
    stopifnot(length(type) == 1)
    stopifnot(is.logical(download) & length(download) == 1)

    ## Use table from the package
    url_table <- recount::recount_url

    ## URLs default to version 2 (disjoint exons).
    if(version == 1) url_table$url <- gsub('v2/', '', url_table$url)

    ## Subset url data
    url_table <- url_table[url_table$project == project, ]
    if(nrow(url_table) == 0) {
        stop("Invalid 'project' argument. There's no such 'project' in the recount_url data.frame.")
    }

    ## If all, download each type individually
    if(type == 'all') {
        urls <- sapply(c(
            'rse-gene', 'rse-exon', 'rse-jx', 'rse-tx',
            'counts-gene', 'counts-exon', 'counts-jx',
            'phenotype', 'files-info', 'samples', 'mean'), function(file_type) {
            Sys.sleep(round(runif(1, 2, 5), 0))
            download_study(project = project, type = file_type,
                outdir = outdir, download = download, version = version, ...)
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
            'rse-jx' = 'rse_jx.Rdata',
            'rse-tx' = 'rse_tx.RData',
            'counts-gene' = 'counts_gene.tsv.gz',
            'counts-exon' = 'counts_exon.tsv.gz',
            'counts-jx' = 'counts_jx.tsv.gz',
            phenotype = paste0(project, '.tsv'),
            'files-info' = 'files_info.tsv',
            mean = paste0('mean_', project, '.bw'),
            'rse-fc' = paste0('rse_fc_', project, '.Rdata')
        )
        url <- url_table$url[url_table$file_name == filename]
        if(length(url) == 0) return(NULL)
        if(download) {
            message(paste(Sys.time(), 'downloading file', filename, 'to',
                outdir))
            xx <- download_retry(
                url = url,
                destfile = file.path(outdir, filename),
                mode = 'wb',
                ...
            )
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
                download_retry(
                    url = sample_urls$url[i],
                    destfile = file.path(outdir, sample_urls$file_name[i]),
                    mode = 'wb',
                    ...
                )
            }, ...)
        }
    }

    ## Return the actual url(s)
    return(invisible(url))
}
