#' Given a set of regions for a chromosome, compute the coverage matrix for a
#' given SRA study.
#'
#' Given a set of genomic regions as created by \link{expressed_regions}, this
#' function computes the coverage matrix for a library size of 40 million 100 bp
#' reads for a given SRA study.
#'
#' @inheritParams expressed_regions
#' @param regions A \link[GenomicRanges]{GRanges-class} object with regions
#' for \code{chr} for which to calculate the coverage matrix.
#' @param chunksize A single integer vector defining the chunksize to use for
#' computing the coverage matrix. Regions will be split into different chunks
#' which can be useful when using a parallel instance as defined by 
#' \code{bpparam}.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} instance which
#' will be used to calculate the coverage matrix in parallel. By default, 
#' \link[BiocParallel]{SerialParam-class} will be used.
#' @param verboseLoad If \code{TRUE} basic status updates for loading the data
#' will be printed.
#' 
#'
#' @return A matrix with one row per region and one column per sample. The
#' numbers in the cells are the counts (number of reads, or fraction in some
#' cases) overlapping the region.
#'
#' @details When using \code{outdir = NULL} the information will be accessed
#' from the web on the fly. If you encounter internet access problems, it might
#' be best to first download the BigWig files using \link{download_study}. This
#' might be the best option if you are accessing all chromosomes for a given
#' project and/or are thinking of using different sets of \code{regions} (for
#' example, from different cutoffs applied to \link{expressed_regions}).
#' Alternatively check the \code{SciServer} section on the vignette to see
#' how to access all the recount data via a R Jupyter Notebook.
#'
#' If you have \code{bwtool} installed, you can use
#' \url{https://github.com/LieberInstitute/recount.bwtool} for faster results.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom utils read.table
#' @import derfinder GenomicRanges RCurl BiocParallel
#'
#' @seealso \link{download_study}, \link[derfinder]{findRegions},
#' \link[derfinder]{railMatrix}
#'
#' @examples
#'
#' if(.Platform$OS.type != 'windows') {
#' ## Reading BigWig files is not supported by rtracklayer on Windows
#'     ## Define expressed regions for study DRP002835, chrY
#'     regions <- expressed_regions('DRP002835', 'chrY', cutoff = 5L, 
#'         maxClusterGap = 3000L)
#'
#'     ## Now calculate the coverage matrix for this study
#'     coverageMatrix <- coverage_matrix('DRP002835', 'chrY', regions)
#'
#'     ## One row per region
#'     identical(length(regions), nrow(coverageMatrix))
#' }
#'

coverage_matrix <- function(project, chr, regions, chunksize = 1000, bpparam = NULL, outdir = NULL, chrlen = NULL, verbose = TRUE, verboseLoad = verbose, ...) {    
    ## Check inputs
    stopifnot(is.character(project) & length(project) == 1)
    stopifnot(is.character(chr) & length(chr) == 1)
    stopifnot((is.numeric(chunksize) | is.integer(chunksize)) & length(chunksize) == 1)
    
    ## Windows-specific info
    if(.Platform$OS.type == 'windows') {
        warning('rtracklayer does not support importing BigWig files on Windows, so this function might not work')
    }
    
    ## Use table from the package
    url_table <- recount::recount_url
    
    ## Subset url data
    url_table <- url_table[url_table$project == project, ]
    stopifnot(nrow(url_table) > 0)
    
    ## Find chromosome length if absent
    if(is.null(chrlen)) {
        chrinfo <- read.table('https://raw.githubusercontent.com/nellore/runs/master/gtex/hg38.sizes',
            col.names = c('chr', 'size'), colClasses = c('character',
            'integer'))
        chrlen <- chrinfo$size[chrinfo$chr == chr]
        stopifnot(length(chrlen) == 1)
    }
    
    samples_i <- which(grepl('[.]bw$', url_table$file_name) & !grepl('mean',
        url_table$file_name))
    ## Check if data is present, otherwise download it
    if(!is.null(outdir)) {
        ## Check sample files
        sampleFiles <- sapply(samples_i, function(i) {
            file.path(outdir, 'bw', url_table$file_name[i])
        })
        if(any(!file.exists(sampleFiles))) {
            download_study(project = project, type = 'samples', outdir = outdir,
                download = TRUE, ...)
        }
        
        ## Check phenotype data
        phenoFile <- file.path(outdir, paste0(project, '.tsv'))
        if(!file.exists(phenoFile)) {
            download_study(project = project, type = 'phenotype',
                outdir = outdir, download = TRUE, ...)
        }
    } else {
        sampleFiles <- download_study(project = project, type = 'samples',
            download = FALSE)
        phenoFile <- download_study(project = project, type = 'phenotype',
            download = FALSE)
    }
        
    ## Read pheno data
    pheno <- .read_pheno(phenoFile, project)
    
    ## Get sample names
    m <- match(url_table$file_name[samples_i], pheno$bigwig_file)
    if(project != 'TCGA') {
        names(sampleFiles) <- pheno$run[m]
    } else {
        names(sampleFiles) <- pheno$gdc_file_id[m]
    }
    
    
    ## Define library size normalization factor
    targetSize <- 40e6
    totalMapped <- pheno$auc[m]
    mappedPerXM <- totalMapped / targetSize
    
    ## Keep only regions from the chr in question
    regions <- regions[seqnames(regions) == chr]
    
    ## Split regions into chunks
    nChunks <- length(regions) %/% chunksize
    if(length(regions) %% chunksize > 0) nChunks <- nChunks + 1
    names(regions) <- seq_len(length(regions))
    
    ## Split regions into chunks
    if(nChunks == 1) {
        regs_split <- list(regions)
        names(regs_split) <- '1'
    } else {
        regs_split <- split(regions, cut(seq_len(length(regions)),
            breaks = nChunks, labels = FALSE))
    }
    
    ## Define bpparam
    if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
    
    ## Load coverage data
    resChunks <- lapply(regs_split, derfinder:::.railMatChrRegion,
        sampleFiles = sampleFiles, chr = chr, mappedPerXM = mappedPerXM,
        L = 1, verbose = verbose, BPPARAM.railChr = bpparam,
        verboseLoad = verboseLoad, chrlen = chrlen)
    
    ## Group results from chunks
    coverageMatrix <- do.call(rbind, lapply(resChunks, '[[', 'coverageMatrix'))
    
    ## Finish
    return(coverageMatrix)
}


## Helper function for reading the phenotype files
.read_pheno <- function(phenoFile, project) {
    if(project %in% c('SRP012682', 'TCGA')) {
        subsets <- c('SRP012682' = 'gtex', 'TCGA' = 'tcga')
        res <- all_metadata(subsets[project], verbose = FALSE)
    } else { 
        info <- readLines(phenoFile)
        res <- read.table(text = info[grepl(paste0('^project|^', project),
            info)], header = TRUE, stringsAsFactors = FALSE, sep = '\t',
            comment.char = '', quote = '')
    }
    return(res)
}
