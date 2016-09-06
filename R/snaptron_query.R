#' Query Snaptron to get data from exon-exon junctions present in Intropolis
#'
#' This function uses the Snaptron API to query specific exon-exon junctions
#' that are available via Intropolis as described in the vignette.
#'
#' @param junctions A \link[GenomicRanges]{GRanges-class} object with the
#' exon-exon junctions of interest. The chromosome names should be in UCSC
#' format, such as 'chr1'. The strand information is ignored in the query.
#' @param version Either \code{1} or \code{2}. Version 1 of Intropolis has the
#' exon-exon junctions from about 20 thousand RNA-seq samples in hg19 
#' coordinates. Version 2 has the data from about 50 thousand RNA-seq samples #' aligned to hg38.
#' @param verbose If \code{TRUE} status updates will be printed.
#'
#' @return A \link[GenomicRanges]{GRanges-class} object with the results from
#' the Snaptron query. For information on the different columns please see
#' \url{http://stingray.cs.jhu.edu:8090/snaptron/docs/}.
#'
#' @references Please cite \url{http://stingray.cs.jhu.edu:8090/snaptron/docs/}
#' if you use this function as Snaptron is a separate project from recount.
#' Thank you!
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @import GenomicRanges
#' @import RCurl
#'
#' @examples
#' 
#' library('GenomicRanges')
#' ## Define some exon-exon junctions (hg19 coordinates)
#' junctions <- GRanges(seqnames = 'chr2', IRanges(
#'     start = c(28971710:28971712, 29555081:29555083, 29754982:29754984),
#'     end = c(29462417:29462419, 29923338:29923340, 29917714:29917716)))
#'
#' ## Check against Snaptron version 1 (hg19 coordinates)
#' snaptron_query(junctions)
#'
#' ## Check another set of junctions against version 2 (more data, hg38 
#'    coordinates)
#' junctions_v2 <- GRanges(seqnames = 'chr2', IRanges(
#'     start = 29532116:29532118, end = 29694848:29694850))
#' snaptron_query(junctions_v2, version = 2) 
#'

snaptron_query <- function(junctions, version = 1, verbose = TRUE) {
    ## Check input
    stopifnot(is(junctions, 'GRanges'))
    stopifnot(all(grepl('chr', seqlevels(junctions))))
    stopifnot(version %in% c(1, 2))
    ver <- ifelse(version == 1, 'srav1', 'srav2')
    
    
    ## Build query URLs
    urls <- paste0('http://stingray.cs.jhu.edu:8090/', ver,
        '/snaptron?regions=', seqnames(junctions), ':', start(junctions), '-',
        end(junctions), '&exact=1&header=0')
    
    ## Get results
    if(verbose) message(paste(Sys.time(), 'querying Snaptron'))
    query_res <- getURL(urls)
    
    ## Split by line
    if(verbose) message(paste(Sys.time(), 'processing results'))
    query_split <- strsplit(query_res, '\n')
    names(query_split) <- NULL
    
    ## Are there any valid ones?
    valid <- which(elementNROWS(query_split) > 0)
    if(length(valid) == 0) {
        message(paste(Sys.time(),
            'found no exon-exon junctions in Intropolis version', version,
            'matching your query: this version uses',
            ifelse(version == 1, 'hg19', 'hg38'), 'coordinates.'))
        return(NULL)
    }
    
    ## Extract information
    if(verbose) message(paste(Sys.time(), 'extracting information'))
    info <- lapply(query_split[valid], function(jxs) { 
        matrix(strsplit(jxs, '\t')[[1]], ncol = 19)
    })
    
    ## Build result
    res <- do.call(rbind, info)
    colnames(res) <- c('type', 'snaptron_id', 'chromosome', 'start', 'end',
        'length', 'strand', 'annotated', 'left_motif', 'right_motif',
        'left_annotated', 'right_annotated', 'samples',
        'read_coverage_by_sample', 'samples_count', 'coverage_sum',
        'coverage_avg', 'coverage_median', 'source_dataset_id')
    
    ## Helper function for some special variables
    to_chr_list <- function(x) {
        r <- strsplit(x, ',')
        i <- which(sapply(r, function(y) { y[[1]] == '0' }))
        if(length(i) > 0) r[i] <- NA
        return(CharacterList(r))
    }
    
    ## Build into a GRanges object
    result <- GRanges(seqnames = res[, 'chromosome'], 
        IRanges(as.numeric(res[, 'start']), as.numeric(res[, 'end'])),
        strand = res[, 'strand'])
    
    ## Add other variables
    result$type <- as.factor(res[, 'type'])
    result$snaptron_id <- as.integer(res[, 'snaptron_id'])
    result$annotated <- to_chr_list(res[, 'annotated'])
    result$left_motif <- res[, 'left_motif']
    result$right_motif <- res[, 'right_motif']
    result$left_annotated <- to_chr_list(res[, 'left_annotated'])
    result$right_annotated <- to_chr_list(res[, 'right_annotated'])
    result$samples <- IntegerList(strsplit(res[, 'samples'], ','))
    result$read_coverage_by_sample <- IntegerList(strsplit(res[,
        'read_coverage_by_sample'], ','))
    result$samples_count <- as.integer(res[, 'samples_count'])
    result$coverage_sum <- as.integer(res[, 'coverage_sum'])
    result$coverage_avg <- as.numeric(res[, 'coverage_avg'])
    result$coverage_median <- as.numeric(res[, 'coverage_median'])
    result$source_dataset_id <- as.integer(res[, 'source_dataset_id'])
    
    ## Finish
    return(result)
}
