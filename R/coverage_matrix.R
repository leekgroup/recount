#' Given a set of regions for a chromosome, compute the coverage matrix for a
#' given SRA study.
#'
#' Given a set of genomic regions as created by [expressed_regions], this
#' function computes the coverage matrix for a library size of 40 million 100 bp
#' reads for a given SRA study.
#'
#' @inheritParams expressed_regions
#' @param regions A [GRanges-class][GenomicRanges::GRanges-class] object with regions
#' for `chr` for which to calculate the coverage matrix.
#' @param chunksize A single integer vector defining the chunksize to use for
#' computing the coverage matrix. Regions will be split into different chunks
#' which can be useful when using a parallel instance as defined by
#' `bpparam`.
#' @param bpparam A [BiocParallelParam-class][BiocParallel::BiocParallelParam-class] instance which
#' will be used to calculate the coverage matrix in parallel. By default,
#' [SerialParam-class][BiocParallel::SerialParam-class] will be used.
#' @param verboseLoad If `TRUE` basic status updates for loading the data
#' will be printed.
#' @param scale If `TRUE`, the coverage counts will be scaled to read
#' counts based on a library size of 40 million reads. Set `scale` to
#' `FALSE` if you want the raw coverage counts. The scaling method is by
#' AUC, as in the default option of [scale_counts].
#' @param round If `TRUE`, the counts are rounded to integers. Set to
#' `TRUE` if you want to match the defaults of [scale_counts].
#'
#'
#' @return A [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#' object with the counts stored in the assays slot.
#'
#' @details When using `outdir = NULL` the information will be accessed
#' from the web on the fly. If you encounter internet access problems, it might
#' be best to first download the BigWig files using [download_study]. This
#' might be the best option if you are accessing all chromosomes for a given
#' project and/or are thinking of using different sets of `regions` (for
#' example, from different cutoffs applied to [expressed_regions]).
#' Alternatively check the `SciServer` section on the vignette to see
#' how to access all the recount data via a R Jupyter Notebook.
#'
#' If you have `bwtool` installed, you can use
#' <https://github.com/LieberInstitute/recount.bwtool> for faster results.
#' Note that you will need to run [scale_counts] after running
#' `coverage_matrix_bwtool()`.
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @importFrom utils read.table
#' @import derfinder GenomicRanges RCurl BiocParallel
#' SummarizedExperiment S4Vectors
#'
#' @seealso [download_study], [findRegions][derfinder::findRegions],
#' [railMatrix][derfinder::railMatrix]
#'
#' @examples
#'
#' if (.Platform$OS.type != "windows") {
#'     ## Reading BigWig files is not supported by rtracklayer on Windows
#'     ## Define expressed regions for study DRP002835, chrY
#'     regions <- expressed_regions("DRP002835", "chrY",
#'         cutoff = 5L,
#'         maxClusterGap = 3000L
#'     )
#'
#'     ## Now calculate the coverage matrix for this study
#'     rse <- coverage_matrix("DRP002835", "chrY", regions)
#'
#'     ## One row per region
#'     identical(length(regions), nrow(rse))
#' }
coverage_matrix <- function(
        project, chr, regions, chunksize = 1000,
        bpparam = NULL, outdir = NULL, chrlen = NULL, verbose = TRUE,
        verboseLoad = verbose, scale = TRUE, round = FALSE, ...) {
    ## Check inputs
    stopifnot(is.character(project) & length(project) == 1)
    stopifnot(is.character(chr) & length(chr) == 1)
    stopifnot((is.numeric(chunksize) | is.integer(chunksize)) & length(chunksize) == 1)

    ## Windows-specific info
    if (.Platform$OS.type == "windows") {
        warning("rtracklayer does not support importing BigWig files on Windows, so this function might not work")
    }

    ## Use table from the package
    url_table <- recount::recount_url

    ## Subset url data
    url_table <- url_table[url_table$project == project, ]
    if (nrow(url_table) == 0) {
        stop("Invalid 'project' argument. There's no such 'project' in the recount_url data.frame.")
    }

    ## Find chromosome length if absent
    if (is.null(chrlen)) {
        chrinfo <- read.table("https://raw.githubusercontent.com/nellore/runs/master/gtex/hg38.sizes",
            col.names = c("chr", "size"), colClasses = c(
                "character",
                "integer"
            )
        )
        chrlen <- chrinfo$size[chrinfo$chr == chr]
        stopifnot(length(chrlen) == 1)
    }

    samples_i <- which(grepl("[.]bw$", url_table$file_name) & !grepl(
        "mean",
        url_table$file_name
    ))
    ## Check if data is present, otherwise download it
    if (!is.null(outdir)) {
        ## Check sample files
        sampleFiles <- sapply(samples_i, function(i) {
            file.path(outdir, "bw", url_table$file_name[i])
        })
        if (any(!file.exists(sampleFiles))) {
            download_study(
                project = project, type = "samples", outdir = outdir,
                download = TRUE, ...
            )
        }

        ## Check phenotype data
        phenoFile <- file.path(outdir, paste0(project, ".tsv"))
        if (!file.exists(phenoFile)) {
            download_study(
                project = project, type = "phenotype",
                outdir = outdir, download = TRUE, ...
            )
        }
    } else {
        sampleFiles <- download_study(
            project = project, type = "samples",
            download = FALSE
        )
        phenoFile <- download_study(
            project = project, type = "phenotype",
            download = FALSE
        )
    }

    ## Read pheno data
    pheno <- .read_pheno(phenoFile, project)

    ## Get sample names
    m <- match(url_table$file_name[samples_i], pheno$bigwig_file)
    ## Match the pheno with the samples
    pheno <- pheno[m, ]
    if (project != "TCGA") {
        names(sampleFiles) <- pheno$run
    } else {
        names(sampleFiles) <- pheno$gdc_file_id
    }

    ## Define library size normalization factor
    if (scale) {
        targetSize <- 40e6
        totalMapped <- pheno$auc[m]
        mappedPerXM <- totalMapped / targetSize
    } else {
        mappedPerXM <- rep(1, nrow(pheno))
    }


    ## Keep only regions from the chr in question
    regions <- regions[seqnames(regions) == chr]

    ## Split regions into chunks
    nChunks <- length(regions) %/% chunksize
    if (length(regions) %% chunksize > 0) nChunks <- nChunks + 1
    names(regions) <- seq_len(length(regions))

    ## Split regions into chunks
    if (nChunks == 1) {
        regs_split <- list(regions)
        names(regs_split) <- "1"
    } else {
        regs_split <- split(regions, cut(seq_len(length(regions)),
            breaks = nChunks, labels = FALSE
        ))
    }

    ## Define bpparam
    if (is.null(bpparam)) bpparam <- BiocParallel::SerialParam()

    ## Load coverage data
    resChunks <- lapply(regs_split, derfinder:::.railMatChrRegion,
        sampleFiles = sampleFiles, chr = chr, mappedPerXM = mappedPerXM,
        L = 1, verbose = verbose, BPPARAM.railChr = bpparam,
        verboseLoad = verboseLoad, chrlen = chrlen
    )

    ## Group results from chunks
    coverageMatrix <- do.call(rbind, lapply(resChunks, "[[", "coverageMatrix"))

    if (round) coverageMatrix <- round(coverageMatrix, 0)

    ## Build a RSE object
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays = list("counts" = coverageMatrix),
        colData = DataFrame(pheno), rowRanges = regions
    )

    ## Finish
    return(rse)
}


## Helper function for reading the phenotype files
.read_pheno <- function(phenoFile, project) {
    if (project %in% c("SRP012682", "TCGA")) {
        subsets <- c("SRP012682" = "gtex", "TCGA" = "tcga")
        res <- all_metadata(subsets[project], verbose = FALSE)
    } else {
        info <- readLines(phenoFile)
        res <- read.table(
            text = info[grepl(
                paste0("^project|^", project),
                info
            )], header = TRUE, stringsAsFactors = FALSE, sep = "\t",
            comment.char = "", quote = ""
        )
    }
    return(res)
}
