## ----docSetup, bootstrap.show.code = FALSE, dev = device, bootstrap.show.message=FALSE----
## knitrBoostrap and device chunk options
load_install('knitr')
opts_chunk$set(bootstrap.show.code = FALSE, dev = device)
if(!outputIsHTML) opts_chunk$set(bootstrap.show.code = FALSE, dev = device, echo = FALSE)

## ----setup, bootstrap.show.message=FALSE---------------------------------
#### Libraries needed

## Bioconductor
load_install('DESeq2')
if(isEdgeR) load_install('edgeR')

## CRAN
load_install('ggplot2')
if(!is.null(theme)) theme_set(theme)
load_install('knitr')
if(is.null(colors)) {
    load_install('RColorBrewer')
}
load_install('pheatmap')
load_install('DT')
load_install('devtools')

## Working behind the scenes
# load_install('knitcitations')
# load_install('rmarkdown')
## Optionally
# load_install('knitrBootstrap')

#### Code setup

## For ggplot
res.df <- as.data.frame(res)

## Sort results by adjusted p-values
ord <- order(res.df$padj, decreasing = FALSE)
res.df <- res.df[ord, ]
res.df <- cbind(data.frame(Feature = rownames(res.df)), res.df)
rownames(res.df) <- NULL

## ----'PCA'---------------------------------------------------------------
## Transform count data
rld <- rlog(dds)

## Perform PCA analysis and make plot
plotPCA(rld, intgroup = intgroup)

## Get percent of variance explained
data_pca <- plotPCA(rld, intgroup = intgroup, returnData = TRUE)
percentVar <- round(100 * attr(data_pca, "percentVar"))

## ----'sampleDist'--------------------------------------------------------
## Obtain the sample euclidean distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
## Add names based on intgroup
rownames(sampleDistMatrix) <- apply(as.data.frame(colData(rld)[, intgroup]), 1,
    paste, collapse = ' : ')
colnames(sampleDistMatrix) <- NULL

## Define colors to use for the heatmap if none were supplied
if(is.null(colors)) {
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
}

## Make the heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists, color = colors)

## ----'MAplotalpha'-------------------------------------------------------
## MA plot with alpha used in DESeq2::results()
plotMA(res, alpha = metadata(res)$alpha, main = paste('MA plot with alpha =',
    metadata(res)$alpha))

## ----'MAplotalphaHalf'---------------------------------------------------
## MA plot with alpha = 1/2 of the alpha used in DESeq2::results()
plotMA(res, alpha = metadata(res)$alpha / 2,
    main = paste('MA plot with alpha =', metadata(res)$alpha / 2))

## ----'MAplotalpha-nBest'-------------------------------------------------
## MA plot with alpha corresponding to the one that gives the nBest features
nBest.actual <- min(nBest, nrow(head(res.df, n = nBest)))
nBest.alpha <- head(res.df, n = nBest)$padj[nBest.actual]
plotMA(res, alpha = nBest.alpha * 1.00000000000001,
    main = paste('MA plot for top', nBest.actual, 'features'))

## ----pvalueHistogram-----------------------------------------------------
## P-value histogram plot
ggplot(res.df[!is.na(res.df$pvalue), ], aes(x = pvalue)) +
    geom_histogram(alpha=.5, position='identity', bins = 50) +
    labs(title='Histogram of unadjusted p-values') +
    xlab('Unadjusted p-values') +
    xlim(c(0, 1.0005))

## ----pvalueSumm----------------------------------------------------------
## P-value distribution summary
summary(res.df$pvalue)

## ----pvalueTable, results = 'asis'---------------------------------------
## Split features by different p-value cutoffs
pval_table <- lapply(c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1), function(x) {
    data.frame('Cut' = x, 'Count' = sum(res.df$pvalue <= x, na.rm = TRUE))
})
pval_table <- do.call(rbind, pval_table)
if(outputIsHTML) {
    kable(pval_table, format = 'markdown', align = c('c', 'c'))
} else {
    kable(pval_table)
}

## ----padjHistogram-------------------------------------------------------
## Adjusted p-values histogram plot
ggplot(res.df[!is.na(res.df$padj), ], aes(x = padj)) +
    geom_histogram(alpha=.5, position='identity', bins = 50) +
    labs(title=paste('Histogram of', elementMetadata(res)$description[grep('adjusted', elementMetadata(res)$description)])) +
    xlab('Adjusted p-values') +
    xlim(c(0, 1.0005))

## ----padjSumm------------------------------------------------------------
## Adjusted p-values distribution summary
summary(res.df$padj)

## ----padjTable, results = 'asis'-----------------------------------------
## Split features by different adjusted p-value cutoffs
padj_table <- lapply(c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1), function(x) {
    data.frame('Cut' = x, 'Count' = sum(res.df$padj <= x, na.rm = TRUE))
})
padj_table <- do.call(rbind, padj_table)
if(outputIsHTML) {
    kable(padj_table, format = 'markdown', align = c('c', 'c'))
} else {
    kable(padj_table)
}

## ----'topFeatures', results = 'asis'-------------------------------------
## Add search url if appropriate
if(!is.null(searchURL) & outputIsHTML) {
    res.df$Feature <- paste0('<a href="', searchURL, res.df$Feature, '">',
        res.df$Feature, '</a>')
}

for(i in which(colnames(res.df) %in% c('pvalue', 'padj'))) res.df[, i] <- format(res.df[, i], scientific = TRUE)

if(outputIsHTML) {
    datatable(head(res.df, n = nBest), options = list(pagingType='full_numbers', pageLength=10, scrollX='100%'), escape = FALSE, rownames = FALSE) %>% formatRound(which(!colnames(res.df) %in% c('pvalue', 'padj', 'Feature')), digits)
} else {
    res.df_top <- head(res.df, n = 20)
    for(i in which(!colnames(res.df) %in% c('pvalue', 'padj', 'Feature'))) res.df_top[, i] <- round(res.df_top[, i], digits)
    kable(res.df_top)
}

## ----'plotCounts'--------------------------------------------------------
plotCounts_gg <- function(i, dds, intgroup) {
    group <- if (length(intgroup) == 1) {
        colData(dds)[[intgroup]]
    } else if (length(intgroup) == 2) {
        lvls <- as.vector(t(outer(levels(colData(dds)[[intgroup[1]]]), 
            levels(colData(dds)[[intgroup[2]]]), function(x, 
                y) paste(x, y, sep = " : "))))
        droplevels(factor(apply(as.data.frame(colData(dds)[, 
            intgroup, drop = FALSE]), 1, paste, collapse = " : "), 
            levels = lvls))
    } else {
        factor(apply(as.data.frame(colData(dds)[, intgroup, drop = FALSE]), 
            1, paste, collapse = " : "))
    }
    data <- plotCounts(dds, gene=i, intgroup=intgroup, returnData = TRUE)
    data <- cbind(data, data.frame('group' = group))
    main <- rownames(dds)[i]

    ggplot(data, aes(x = group, y = count)) + geom_point() + ylab('Normalized count') + ggtitle(main) + coord_trans(y = "log10")
}
for(i in head(ord, nBestFeatures)) {
    print(plotCounts_gg(i, dds = dds, intgroup = intgroup))
}

## ----'edgeR-BCV', eval = isEdgeR, echo = isEdgeR & outputIsHTML----------
#  ## Make BCV plot for edgeR results
#  plotBCV(dge)

## ----'edgeR-MDS', eval = isEdgeR, echo = isEdgeR & outputIsHTML----------
#  ## Make MDS plot for edgeR results
#  plotMDS(dge)

## ----child = customCode, eval = hasCustomCode----------------------------
#  NA

## ----thecall, echo=FALSE-------------------------------------------------
theCall

## ----reproducibility1, echo=FALSE----------------------------------------
## Date the report was generated
Sys.time()

## ----reproducibility2, echo=FALSE----------------------------------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)

## ----reproducibility3, echo=FALSE-------------------------------------------------------------------------------------
## Session info
options(width = 120)
session_info()

## ----bibliography, results='asis', echo=FALSE, warning = FALSE--------------------------------------------------------
## Print bibliography
bibliography()

