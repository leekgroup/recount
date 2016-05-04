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
    BiocStyle = citation('BiocStyle'),
    derfinder = citation('derfinder')[1], 
    DESeq2 = citation('DESeq2'),
    devtools = citation('devtools'),
    GenomicFeatures = citation('GenomicFeatures'),
    GenomicRanges = citation('GenomicRanges'),
    IRanges = citation('IRanges'),
    knitcitations = citation('knitcitations'),
    knitr = citation('knitr')[3],
    R = citation(),
    regionReport = citation('regionReport'),
    rmarkdown = citation('rmarkdown'),
    rtracklayer = citation('rtracklayer'),
    SummarizedExperiment = citation('SummarizedExperiment'),
    testthat = citation('testthat'),
    TxDb.Hsapiens.UCSC.hg38.knownGene = citation('TxDb.Hsapiens.UCSC.hg38.knownGene')
)

write.bibtex(bibs,
    file = 'quickstartRef.bib')
bib <- read.bibtex('quickstartRef.bib')

## Assign short names
names(bib) <- names(bibs)

## Working on Windows?
windowsFlag <- .Platform$OS.type == 'windows'

## ----'ultraQuick', eval = FALSE------------------------------------------
#  ## Load libraries
#  library('recount')
#  library('SummarizedExperiment')
#  
#  ## Find a project of interest
#  project_info <- abstract_search('GSE32465')
#  
#  ## Browse the project at SRA
#  browse_study(project_info$project)
#  
#  ## Define groups based on the project information available at
#  ## http://www.ncbi.nlm.nih.gov/sra/?term=SRP009615 (found via the "Experiments"
#  ## link) from the SRA website.
#  sample_info <- data.frame(
#      accession = c('SRX110461', 'SRX110462', 'SRX110463', 'SRX110464',
#          'SRX111299', 'SRX111300', 'SRX111301', 'SRX111302', 'SRX111303',
#          'SRX111304', 'SRX111305', 'SRX111306'),
#      group = rep(c('uninduced', 'induced'), 6),
#      gene_target = rep(c('SRF', 'EGR1', 'ATF3'), each = 4),
#      replicate = factor(rep(rep(c(1, 2), each = 2), 3))
#  )
#  
#  ## Download the gene-level RangedSummarizedExperiment data
#  download_study(project_info$project)
#  
#  ## Load the data
#  load(file.path(project_info$project), 'rse_gene.Rdata')
#  
#  ## Scale counts by taking into account the total coverage per sample
#  rse <- scale_counts(rse_gene)
#  
#  ## Add sample information for DE analysis
#  i <- match(colData(rse)$experiment, sample_info$accession)
#  colData(rse)$group <- sample_info$group[i]
#  colData(rse)$gene_target <- sample_info$gene_target[i]
#  colData(rse)$replicate <- sample_info$replicate[i]
#  
#  ## Perform differential expression analysis with DESeq2
#  library('DESeq2')
#  
#  ## Specify design and switch to DESeq2 format
#  dds <- DESeqDataSet(rse, ~ group + gene_target + replicate)
#  
#  ## Perform DE analysis
#  dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target + replicate,
#      fitType = 'local')
#  res <- results(dds)
#  
#  ## Explore results
#  plotMA(res, main="DESeq2 results for SRP009615")
#  
#  ## Make a report with the results
#  library('regionReport')
#  DESeq2Report(dds, res = res, project = 'SRP009615',
#      intgroup = c('group', 'gene_target', 'replicate'), outdir = '.',
#      output = 'SRP009615-results')

## ----'start', message=FALSE----------------------------------------------
## Load libraries
library('recount')
library('SummarizedExperiment')

## ----'search_abstract'---------------------------------------------------
## Find a project of interest
project_info <- abstract_search('GSE32465')

## Explore info
project_info

## ----'browse'------------------------------------------------------------
## Browse the project at SRA
browse_study(project_info$project)

## ----'sample_info'-------------------------------------------------------
## Define groups based on the project information available at
## http://www.ncbi.nlm.nih.gov/sra/?term=SRP009615 (found via the "Experiments"
## link) from the SRA website.
sample_info <- data.frame(
    accession = c('SRX110461', 'SRX110462', 'SRX110463', 'SRX110464',
        'SRX111299', 'SRX111300', 'SRX111301', 'SRX111302', 'SRX111303',
        'SRX111304', 'SRX111305', 'SRX111306'),
    group = rep(c('uninduced', 'induced'), 6),
    gene_target = rep(c('SRF', 'EGR1', 'ATF3'), each = 4),
    replicate = factor(rep(rep(c(1, 2), each = 2), 3))
)

## ----'download', eval = FALSE--------------------------------------------
#  ## Download the gene-level RangedSummarizedExperiment data
#  download_study(project_info$project)
#  
#  ## Load the data
#  load(file.path(project_info$project), 'rse_gene.Rdata')

## ----'download-hidden', echo = FALSE-------------------------------------
## Currently the data is not uploaded yet
rse_gene <- rse_gene_SRP009615

## ----'explore_rse'-------------------------------------------------------
rse_gene

## This is the sample phenotype data provided by the recount project
colData(rse_gene)

## At the gene level, the row data includes the names of the genes and
## the sum of the reduced exons widths, which can be used for taking into
## account the gene length.
rowData(rse_gene)

## ----'scale_counts'------------------------------------------------------
## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)

## ----'add_sample_info'---------------------------------------------------
## Add sample information for DE analysis
i <- match(colData(rse)$experiment, sample_info$accession)
colData(rse)$group <- sample_info$group[i]
colData(rse)$gene_target <- sample_info$gene_target[i]
colData(rse)$replicate <- sample_info$replicate[i]

## ----'de_analysis'-------------------------------------------------------
## Perform differential expression analysis with DESeq2
library('DESeq2')

## Specify design and switch to DESeq2 format
dds <- DESeqDataSet(rse, ~ group + gene_target + replicate)

## Perform DE analysis
dds <- DESeq(dds, test = 'LRT', reduced = ~ gene_target + replicate,
    fitType = 'local')
res <- results(dds)

## ----'ma_plot'-----------------------------------------------------------
## Explore results
plotMA(res, main="DESeq2 results for SRP009615")

## ----'make_report', eval = FALSE-----------------------------------------
#  ## Make a report with the results
#  library('regionReport')
#  report <- DESeq2Report(dds, res = res, project = 'SRP009615',
#      intgroup = c('group', 'gene_target', 'replicate'), outdir = '.',
#      output = 'SRP009615-results')

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
    intgroup = c('group', 'gene_target', 'replicate'), outdir = '.',
    output = 'SRP009615-results', device = 'png', template = 'SRP009615-results-template.Rmd')
    
## Clean up
file.remove('SRP009615-results-template.Rmd')

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
#  load('../../runs/recount2/genes/ucsc-knowngene-hg38-exons.Rdata')
#  recount_exons <- exons
#  save(recount_exons, file = '../data/recount_exons.RData')
#  load('../.../runs/recount2/genes/ucsc-knowngene-hg38-genes-bp-length.Rdata')
#  recount_genes <- genes
#  save(recount_genes, file = '../data/recount_genes.RData')
#  
#  ## URL table
#  load('../../runs/recount2/fileinfo/upload_table.Rdata')
#  recount_url <- upload_table
#  ## Fake urls for now
#  recount_url$url <- 'https://lcolladotor.shinyapps.io/recount/ucsc-knowngene-hg38-genes-bp-length.Rdata'
#  save(recount_url, file = '../data/recount_url.RData')
#  
#  ## Abstract info
#  load('../../runs/recount2/metadata_web/meta_web_sra.Rdata')
#  recount_abstract <- meta_web[, 2:4]
#  recount_abstract$project <- gsub('.*">|</a>', '', meta_web$accession)
#  save(recount_abstract, file = '../data/recount_abstract.RData')
#  
#  ## Example rse_gene file
#  system('scp e:/dcl01/leek/data/gtex_work/runs/recount2/rse/rse_sra/SRP009615/rse_gene.Rdata .')
#  load('rse_gene.Rdata')
#  rse_gene_SRP009615 <- rse_gene
#  save(rse_gene_SRP009615, file = '../data/rse_gene_SRP009615.RData')
#  unlink('rse_gene.Rdata')

## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE--------------------------------------------------
## Print bibliography
bibliography()

