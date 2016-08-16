<a href="http://www.bioconductor.org/packages/release/bioc/html/recount.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/recount.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="https://bioconductor.org/packages/stats/bioc/recount/"><img border="0" src="http://www.bioconductor.org/shields/downloads/recount.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/recount/"><img border="0" src="http://www.bioconductor.org/shields/posts/recount.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/release/bioc/html/recount.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/recount.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

Status: Travis CI [![Build Status](https://travis-ci.org/leekgroup/recount.svg?branch=master)](https://travis-ci.org/leekgroup/recount),
Bioc-release <a href="http://www.bioconductor.org/packages/release/bioc/html/recount.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/release/recount.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/recount/"><img border="0" src="http://www.bioconductor.org/shields/build/release/bioc/recount.svg" title="build results; click for full report"></a>,
Bioc-devel <a href="http://www.bioconductor.org/packages/devel/bioc/html/recount.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/devel/recount.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/recount/"><img border="0" src="http://www.bioconductor.org/shields/build/devel/bioc/recount.svg" title="build results; click for full report"></a>.

Bioc-release <a href="https://bioconductor.org/developers/how-to/unitTesting-guidelines/#coverage"><img border="0" src="http://www.bioconductor.org/shields/coverage/release/recount.svg" title="Test coverage percentage, or 'unknown'"></a>, Bioc-devel <a href="https://codecov.io/github/Bioconductor-mirror/recount?branch=master"><img border="0" src="http://www.bioconductor.org/shields/coverage/devel/recount.svg" title="Test coverage percentage, or 'unknown'"></a>, Codecov [![codecov.io](https://codecov.io/github/leekgroup/recount/coverage.svg?branch=master)](https://codecov.io/github/leekgroup/recount?branch=master)

recount
=======

Explore and download data from the recount project available at the [recount website](https://jhubiostatistics.shinyapps.io/recount/). Using the `recount` package you can download _RangedSummarizedExperiment_ objects at the gene, exon or exon-exon junctions level, the raw counts, the phenotype metadata used, the urls to the sample coverage bigWig files or the mean coverage bigWig file for a particular study. The _RangedSummarizedExperiment_ objects can be used by different packages for performing differential expression analysis. Using [derfinder](http://bioconductor.org/packages/derfinder) you can perform annotation-agnostic differential expression analyses with the data from the recount project. 

For more information about `recount` check the vignettes.

# Installation instructions

Get R 3.3.0 from [CRAN](http://cran.r-project.org/).

```R
## From Bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite('recount')

## Suggested:
biocLite(c('derfinder', 'DESeq2'))
```

# Vignettes

The vignettes for this package can be viewed via [Bioconductor's website](http://www.bioconductor.org/packages/recount) (manual [backup](http://leekgroup.github.io/recount/)).


# Citation

Below is the citation output from using `citation('recount')` in R. Please 
run this yourself to check for any updates on how to cite __recount__.

To cite the __recount__ package in publications use:

Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD, Jaffe AE, Langmead B and Leek JT (2016). “recount: A large-scale resource of analysis-ready RNA-seq expression data.” _bioRxiv_. doi: 10.1101/068478 (URL: http://doi.org/10.1101/068478), <URL:
http://biorxiv.org/content/early/2016/08/08/068478>.

A BibTeX entry for LaTeX users is

@Article{,
    title = {recount: A large-scale resource of analysis-ready RNA-seq expression data},
    author = {Leonardo Collado-Torres and Abhinav Nellore and Kai Kammers and Shannon E. Ellis and Margaret A. Taub and Kasper D. Hansen and  and Andrew E. Jaffe and Ben Langmead and Jeffrey T. Leek},
    year = {2016},
    journal = {bioRxiv},
	doi = {10.1101/068478}
    url = {http://biorxiv.org/content/early/2016/08/08/068478},
}

# Testing

Testing on Bioc-devel is feasible thanks to [R Travis](http://docs.travis-ci.com/user/languages/r/) as well as Bioconductor's nightly build.
