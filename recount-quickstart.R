## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE-----------
## Track time spent on making the vignette
startTime <- Sys.time()

## Bib setup
library('knitcitations')

## Load knitcitations with a clean bibliography
cleanbib()
cite_options(hyperlink = 'to.doc', citation_format = 'text', style = 'html')
# Note links won't show for now due to the following issue
# https://github.com/cboettig/knitcitations/issues/63

## Write bibliography information
bibs <- c(
    AnnotationDbi = citation('AnnotationDbi'),
    BiocParallel = citation('BiocParallel'),
    BiocStyle = citation('BiocStyle'),
    derfinder = citation('derfinder')[1], 
    DESeq2 = citation('DESeq2'),
    devtools = citation('devtools'),
    downloader = citation('downloader'),
    GEOquery = citation('GEOquery'),
    GenomeInfoDb = citation('GenomeInfoDb'),
    GenomicFeatures = citation('GenomicFeatures'),
    GenomicRanges = citation('GenomicRanges'),
    IRanges = citation('IRanges'),
    knitcitations = citation('knitcitations'),
    knitr = citation('knitr')[3],
    org.Hs.eg.db = citation('org.Hs.eg.db'),
    R = citation(),
    RCurl = citation('RCurl'),
	recount = citation('recount'),
    regionReport = citation('regionReport'),
    rentrez = citation('rentrez'),
    rmarkdown = citation('rmarkdown'),
    rtracklayer = citation('rtracklayer'),
    S4Vectors = citation('S4Vectors'),
    SummarizedExperiment = citation('SummarizedExperiment'),
    testthat = citation('testthat'),
    TxDb.Hsapiens.UCSC.hg38.knownGene = citation('TxDb.Hsapiens.UCSC.hg38.knownGene')
)

write.bibtex(bibs,
    file = 'quickstartRef.bib')
bib <- read.bibtex('quickstartRef.bib')

## Assign short names
names(bib) <- names(bibs)

## ----'ultraQuick', eval = FALSE------------------------------------------
#  ## Load library
#  library('recount')
#  
#  ## Find a project of interest
#  project_info <- abstract_search('GSE32465')
#  
#  ## Download the gene-level RangedSummarizedExperiment data
#  download_study(project_info$project)
#  
#  ## Load the data
#  load(file.path(project_info$project, 'rse_gene.Rdata'))
#  
#  ## Browse the project at SRA
#  browse_study(project_info$project)
#  
#  ## View GEO ids
#  colData(rse_gene)$geo_accession
#  
#  ## Extract the sample characteristics
#  geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)
#  
#  ## Note that the information for this study is a little inconsistent, so we
#  ## have to fix it.
#  geochar <- do.call(rbind, lapply(geochar, function(x) {
#      if('cells' %in% colnames(x)) {
#          colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
#          return(x)
#      } else {
#          return(x)
#      }
#  }))
#  
#  ## We can now define some sample information to use
#  sample_info <- data.frame(
#      run = colData(rse_gene)$run,
#      group = ifelse(grepl('uninduced', colData(rse_gene)$title), 'uninduced', 'induced'),
#      gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
#          'targeting ')[[1]][2], ' gene')[[1]][1] }),
#      cell.line = geochar$cell.line
#  )
#  
#  ## Scale counts by taking into account the total coverage per sample
#  rse <- scale_counts(rse_gene)
#  
#  ## Add sample information for DE analysis
#  colData(rse)$group <- sample_info$group
#  colData(rse)$gene_target <- sample_info$gene_target
#  
#  ## Perform differential expression analysis with DESeq2
#  library('DESeq2')
#  
#  ## Specify design and switch to DESeq2 format
#  dds <- DESeqDataSet(rse, ~ gene_target + group)
#  
#  ## Perform DE analysis
#  dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target, fitType = 'local')
#  res <- results(dds)
#  
#  ## Explore results
#  plotMA(res, main="DESeq2 results for SRP009615")
#  
#  ## Make a report with the results
#  library('regionReport')
#  DESeq2Report(dds, res = res, project = 'SRP009615',
#      intgroup = c('group', 'gene_target'), outdir = '.',
#      output = 'SRP009615-results')

## ----'er_analysis', eval = FALSE-----------------------------------------
#  ## Define expressed regions for study SRP009615, only for chromosome Y
#  regions <- expressed_regions('SRP009615', 'chrY', cutoff = 5L,
#      maxClusterGap = 3000L)
#  
#  ## Compute coverage matrix for study SRP009615, only for chromosome Y
#  system.time( coverageMatrix <- coverage_matrix('SRP009615', 'chrY', regions) )
#  
#  ## Round the coverage matrix to integers
#  covMat <- round(coverageMatrix, 0)
#  
#  ## Get phenotype data for study SRP009615
#  pheno_url <- download_study(project = project_info$project, type = 'phenotype',
#      download = FALSE)
#  pheno <- read.table(pheno_url, header = TRUE, stringsAsFactors = FALSE)
#  
#  ## We can sort the table to make sure everything is in the correct order
#  pheno <- pheno[match(colnames(coverageMatrix), pheno$run), ]
#  
#  ## Complete the phenotype table with the data we got from GEO
#  m <- match(pheno$run, sample_info$run)
#  pheno <- cbind(pheno, sample_info[m, 2:3])
#  
#  ## Build a DESeqDataSet
#  dds_ers <- DESeqDataSetFromMatrix(countData = covMat, colData = pheno,
#      design =  ~ gene_target + group)
#  
#  ## Perform differential expression analysis with DESeq2 at the ER-level
#  dds_ers <- DESeq(dds_ers, test = 'LRT', reduced = ~ gene_target,
#      fitType = 'local')
#  res_ers <- results(dds_ers)
#  
#  ## Explore results
#  plotMA(res_ers, main="DESeq2 results for SRP009615 (ER-level, chrY)")
#  
#  ## Create a more extensive exploratory report
#  DESeq2Report(dds_ers, res = res_ers,
#      project = 'SRP009615 (ER-level, chrY)',
#      intgroup = c('group', 'gene_target'), outdir = '.',
#      output = 'SRP009615-results-ER-level-chrY')

## ----'install', eval = FALSE---------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("recount")

## ----'start', message=FALSE----------------------------------------------
## Load recount R package
library('recount')

## ----'search_abstract'---------------------------------------------------
## Find a project of interest
project_info <- abstract_search('GSE32465')

## Explore info
project_info

## ----'download'----------------------------------------------------------
## Download the gene-level RangedSummarizedExperiment data
download_study(project_info$project)

## Load the data
load(file.path(project_info$project, 'rse_gene.Rdata'))

## ----'explore_rse'-------------------------------------------------------
rse_gene

## This is the sample phenotype data provided by the recount project
colData(rse_gene)

## At the gene level, the row data includes the gene ENTREZ ids, the gene
## symbols and the sum of the reduced exons widths, which can be used for 
## taking into account the gene length.
rowData(rse_gene)

## At the exon level, you can get the gene ENTREZ ids from the names of:
# rowRanges(rse_exon)

## ----'browse'------------------------------------------------------------
## Browse the project at SRA
browse_study(project_info$project)

## ----'sample_info', warning = FALSE--------------------------------------
## View GEO ids
colData(rse_gene)$geo_accession

## Extract the sample characteristics
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)

## Note that the information for this study is a little inconsistent, so we
## have to fix it.
geochar <- do.call(rbind, lapply(geochar, function(x) {
    if('cells' %in% colnames(x)) {
        colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
        return(x)
    } else {
        return(x)
    }
}))

## We can now define some sample information to use
sample_info <- data.frame(
    run = colData(rse_gene)$run,
    group = ifelse(grepl('uninduced', colData(rse_gene)$title), 'uninduced', 'induced'),
    gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
        'targeting ')[[1]][2], ' gene')[[1]][1] }),
    cell.line = geochar$cell.line
)

## ----'scale_counts'------------------------------------------------------
## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)

## ----'add_sample_info'---------------------------------------------------
## Add sample information for DE analysis
colData(rse)$group <- sample_info$group
colData(rse)$gene_target <- sample_info$gene_target

## ----'de_analysis'-------------------------------------------------------
## Perform differential expression analysis with DESeq2
library('DESeq2')

## Specify design and switch to DESeq2 format
dds <- DESeqDataSet(rse, ~ gene_target + group)

## Perform DE analysis
dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target, fitType = 'local')
res <- results(dds)

## ----'ma_plot'-----------------------------------------------------------
## Explore results
plotMA(res, main="DESeq2 results for SRP009615")

## ----'make_report', eval = FALSE-----------------------------------------
#  ## Make a report with the results
#  library('regionReport')
#  report <- DESeq2Report(dds, res = res, project = 'SRP009615',
#      intgroup = c('group', 'gene_target'), outdir = '.',
#      output = 'SRP009615-results', nBest = 10, nBestFeatures = 2)

## ----'make_report_real', echo = FALSE, results = 'hide'------------------
library('regionReport')

## Make it so that the report will be available as a vignette
original <- readLines(system.file('DESeq2Exploration', 'DESeq2Exploration.Rmd',
    package = 'regionReport'))
vignetteInfo <- c(
    'vignette: >',
    '  %\\VignetteEngine{knitr::rmarkdown}',
    '  %\\VignetteIndexEntry{Basic DESeq2 results exploration}',
    '  %\\VignetteEncoding{UTF-8}'
)
new <- c(original[1:12], vignetteInfo, original[13:length(original)])
writeLines(new, 'SRP009615-results-template.Rmd')

## Now create the report
report <- DESeq2Report(dds, res = res, project = 'SRP009615',
    intgroup = c('group', 'gene_target'), outdir = '.',
    output = 'SRP009615-results', device = 'png', template = 'SRP009615-results-template.Rmd', nBest = 10, nBestFeatures = 2)
    
## Clean up
file.remove('SRP009615-results-template.Rmd')

## ----'geneSymbols'-------------------------------------------------------
## Load required library
library('org.Hs.eg.db')

## Extract ENTREZ gene ids
entrez <- names(recount_genes)

## Find the gene information we are interested in
gene_info <- select(org.Hs.eg.db, entrez, c('ENTREZID', 'GENENAME', 'SYMBOL'),
    'ENTREZID')

## Explore part of the results
dim(gene_info)
head(gene_info)

## ----'define_ers', eval = .Platform$OS.type != 'windows'-----------------
## Define expressed regions for study SRP009615, only for chromosome Y
regions <- expressed_regions('SRP009615', 'chrY', cutoff = 5L, 
    maxClusterGap = 3000L)

## Briefly explore the resulting regions
regions

## ----'compute_covMat', eval = .Platform$OS.type != 'windows'-------------
## Compute coverage matrix for study SRP009615, only for chromosome Y
system.time( coverageMatrix <- coverage_matrix('SRP009615', 'chrY', regions) )

## Explore the matrix a bit
dim(coverageMatrix)
head(coverageMatrix)

## ----'to_integer', eval = .Platform$OS.type != 'windows'-----------------
## Round the coverage matrix to integers
covMat <- round(coverageMatrix, 0)

## ----'phenoData', eval = .Platform$OS.type != 'windows'------------------
## Get phenotype data for study SRP009615
pheno_url <- download_study(project = project_info$project, type = 'phenotype',
    download = FALSE)
pheno <- read.table(pheno_url, header = TRUE, stringsAsFactors = FALSE, sep = '\t')

## We can sort the table to make sure everything is in the correct order
pheno <- pheno[match(colnames(coverageMatrix), pheno$run), ]

## Complete the phenotype table with the data we got from GEO
m <- match(pheno$run, sample_info$run)
pheno <- cbind(pheno, sample_info[m, 2:3])

## Explore the phenotype data a little bit
head(pheno)

## ----'ers_dds', eval = .Platform$OS.type != 'windows'--------------------
## Build a DESeqDataSet
dds_ers <- DESeqDataSetFromMatrix(countData = covMat, colData = pheno,
    design =  ~ gene_target + group)

## ----'de_analysis_ers', eval = .Platform$OS.type != 'windows'------------
## Perform differential expression analysis with DESeq2 at the ER-level
dds_ers <- DESeq(dds_ers, test = 'LRT', reduced = ~ gene_target,
    fitType = 'local')
res_ers <- results(dds_ers)

## ----'ma_plot_ers', eval = .Platform$OS.type != 'windows'----------------
## Explore results
plotMA(res_ers, main="DESeq2 results for SRP009615 (ER-level, chrY)")

## ----'report2', eval = FALSE---------------------------------------------
#  ## Create the report
#  report2 <- DESeq2Report(dds_ers, res = res_ers,
#      project = 'SRP009615 (ER-level, chrY)',
#      intgroup = c('group', 'gene_target'), outdir = '.',
#      output = 'SRP009615-results-ER-level-chrY')

## ----snaptron------------------------------------------------------------
library('GenomicRanges')
junctions <- GRanges(seqnames = 'chr2', IRanges(
    start = c(28971711, 29555082, 29754983),
    end = c(29462418, 29923339, 29917715)))

snaptron_query(junctions)

## ----'installDer', eval = FALSE------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("recount")

## ----'citation'----------------------------------------------------------
## Citation info
citation('recount')

## ----createVignette, eval=FALSE------------------------------------------
#  ## Create the vignette
#  library('rmarkdown')
#  system.time(render('recount-quickstart.Rmd', 'BiocStyle::html_document'))
#  
#  ## Extract the R code
#  library('knitr')
#  knit('recount-quickstart.Rmd', tangle = TRUE)

## ----createVignette2-----------------------------------------------------
## Clean up
file.remove('quickstartRef.bib')

## ----reproduce1, echo=FALSE----------------------------------------------
## Date the vignette was generated
Sys.time()

## ----reproduce2, echo=FALSE----------------------------------------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)

## ----reproduce3, echo=FALSE-------------------------------------------------------------------------------------------
## Session info
library('devtools')
options(width = 120)
session_info()

## ----'datasetup', echo = FALSE, eval = FALSE--------------------------------------------------------------------------
#  ## Code for re-creating the data distributed in this package
#  
#  ## Genes/exons
#  library('GenomicRanges')
#  load('../../recount-website/genes/ucsc-knowngene-hg38-exons.Rdata')
#  recount_exons <- exons
#  save(recount_exons, file = '../data/recount_exons.RData')
#  load('../../recount-website/genes/ucsc-knowngene-hg38-genes-bp-length.Rdata')
#  recount_genes <- genes
#  save(recount_genes, file = '../data/recount_genes.RData', compress = 'xz')
#  
#  ## URL table
#  load('../../recount-website/fileinfo/upload_table.Rdata')
#  recount_url <- upload_table
#  ## Fake urls for now
#  is.bw <- grepl('[.]bw$', recount_url$file_name)
#  recount_url$url <- NA
#  recount_url$url[!is.bw] <- paste0('http://duffel.rail.bio/recount/',
#      recount_url$project[!is.bw], '/', recount_url$file_name[!is.bw])
#  recount_url$url[is.bw] <- paste0('http://duffel.rail.bio/recount/',
#      recount_url$project[is.bw], '/bw/', recount_url$file_name[is.bw])
#  save(recount_url, file = '../data/recount_url.RData', compress = 'xz')
#  
#  ## Abstract info
#  load('../../recount-website/website/meta_web.Rdata')
#  recount_abstract <- meta_web[, 2:4]
#  recount_abstract$project <- gsub('.*">|</a>', '', meta_web$accession)
#  Encoding(recount_abstract$abstract) <- 'latin1'
#  recount_abstract$abstract <- iconv(recount_abstract$abstract, 'latin1', 'UTF-8')
#  save(recount_abstract, file = '../data/recount_abstract.RData',
#      compress = 'bzip2')
#  
#  ## Example rse_gene file
#  system('scp e:/dcl01/leek/data/gtex_work/runs/recount2/rse/rse_sra/SRP009615/rse_gene.Rdata .')
#  load('rse_gene.Rdata')
#  rse_gene_SRP009615 <- rse_gene
#  save(rse_gene_SRP009615, file = '../data/rse_gene_SRP009615.RData',
#      compress = 'xz')
#  unlink('rse_gene.Rdata')

## ----'cleanup', echo = FALSE------------------------------------------------------------------------------------------
## Clean up
unlink('SRP009615', recursive = TRUE)

## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE--------------------------------------------------
## Print bibliography
bibliography()

