
<!-- README.md is generated from README.Rmd. Please edit that file -->

# recount <img src="man/figures/logo.png" align="right" width="400px" />

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/recount.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/recount)
[![Codecov test
coverage](https://codecov.io/gh/leekgroup/recount/branch/master/graph/badge.svg)](https://codecov.io/gh/leekgroup/recount?branch=master)
[![R build
status](https://github.com/leekgroup/recount/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/leekgroup/recount/actions)
<!-- badges: end -->

Explore and download data from the recount project available at the
[recount2 website](https://jhubiostatistics.shinyapps.io/recount/).
Using the `recount` package you can download
*RangedSummarizedExperiment* objects at the gene, exon or exon-exon
junctions level, the raw counts, the phenotype metadata used, the urls
to the sample coverage bigWig files or the mean coverage bigWig file for
a particular study. The *RangedSummarizedExperiment* objects can be used
by different packages for performing differential expression analysis.
Using [derfinder](http://bioconductor.org/packages/derfinder) you can
perform annotation-agnostic differential expression analyses with the
data from the recount project.

## Documentation

For more information about `recount` check the vignettes [through
Bioconductor](http://bioconductor.org/packages/recount) or at the
[documentation website](http://leekgroup.github.io/recount).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `recount` using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("recount")
```

## Citation

Below is the citation output from using `citation('recount')` in R.
Please run this yourself to check for any updates on how to cite
**recount**.

``` r
print(citation('recount'), bibtex = TRUE)
#> 
#> Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD,
#> Jaffe AE, Langmead B, Leek JT (2017). "Reproducible RNA-seq analysis
#> using recount2." _Nature Biotechnology_. doi: 10.1038/nbt.3838 (URL:
#> https://doi.org/10.1038/nbt.3838), <URL:
#> http://www.nature.com/nbt/journal/v35/n4/full/nbt.3838.html>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Reproducible RNA-seq analysis using recount2},
#>     author = {Leonardo Collado-Torres and Abhinav Nellore and Kai Kammers and Shannon E. Ellis and Margaret A. Taub and Kasper D. Hansen and Andrew E. Jaffe and Ben Langmead and Jeffrey T. Leek},
#>     year = {2017},
#>     journal = {Nature Biotechnology},
#>     doi = {10.1038/nbt.3838},
#>     url = {http://www.nature.com/nbt/journal/v35/n4/full/nbt.3838.html},
#>   }
#> 
#> Collado-Torres L, Nellore A, Jaffe AE (2017). "recount workflow:
#> Accessing over 70,000 human RNA-seq samples with Bioconductor [version
#> 1; referees: 1 approved, 2 approved with reservations]."
#> _F1000Research_. doi: 10.12688/f1000research.12223.1 (URL:
#> https://doi.org/10.12688/f1000research.12223.1), <URL:
#> https://f1000research.com/articles/6-1558/v1>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {recount workflow: Accessing over 70,000 human RNA-seq samples with Bioconductor [version 1; referees: 1 approved, 2 approved with reservations]},
#>     author = {Leonardo Collado-Torres and Abhinav Nellore and Andrew E. Jaffe},
#>     year = {2017},
#>     journal = {F1000Research},
#>     doi = {10.12688/f1000research.12223.1},
#>     url = {https://f1000research.com/articles/6-1558/v1},
#>   }
#> 
#> Ellis SE, Collado-Torres L, Jaffe AE, Leek JT (2018). "Improving the
#> value of public RNA-seq expression data by phenotype prediction."
#> _Nucl. Acids Res._. doi: 10.1093/nar/gky102 (URL:
#> https://doi.org/10.1093/nar/gky102), <URL:
#> https://doi.org/10.1093/nar/gky102>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Improving the value of public RNA-seq expression data by phenotype prediction},
#>     author = {Shannon E. Ellis and Leonardo Collado-Torres and Andrew E. Jaffe and Jeffrey T. Leek},
#>     year = {2018},
#>     journal = {Nucl. Acids Res.},
#>     doi = {10.1093/nar/gky102},
#>     url = {https://doi.org/10.1093/nar/gky102},
#>   }
#> 
#> Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD,
#> Jaffe AE, Langmead B, Leek JT (2020). _Explore and download data from
#> the recount project_. doi: 10.18129/B9.bioc.recount (URL:
#> https://doi.org/10.18129/B9.bioc.recount),
#> https://github.com/leekgroup/recount - R package version 1.13.2, <URL:
#> http://www.bioconductor.org/packages/recount>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Explore and download data from the recount project},
#>     author = {Leonardo Collado-Torres and Abhinav Nellore and Kai Kammers and Shannon E. Ellis and Margaret A. Taub and Kasper D. Hansen and Andrew E. Jaffe and Ben Langmead and Jeffrey T. Leek},
#>     year = {2020},
#>     url = {http://www.bioconductor.org/packages/recount},
#>     note = {https://github.com/leekgroup/recount - R package version 1.13.2},
#>     doi = {10.18129/B9.bioc.recount},
#>   }
#> 
#> Frazee AC, Langmead B, Leek JT (2011). "ReCount: A multi-experiment
#> resource of analysis-ready RNA-seq gene count datasets." _BMC
#> Bioinformatics_. doi: 10.1186/1471-2105-12-449 (URL:
#> https://doi.org/10.1186/1471-2105-12-449), <URL:
#> https://doi.org/10.1186/1471-2105-12-449>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {ReCount: A multi-experiment resource of analysis-ready RNA-seq gene count datasets},
#>     author = {Alyssa C. Frazee and Ben Langmead and Jeffrey T. Leek},
#>     year = {2011},
#>     journal = {BMC Bioinformatics},
#>     doi = {10.1186/1471-2105-12-449},
#>     url = {https://doi.org/10.1186/1471-2105-12-449},
#>   }
#> 
#> Razmara A, Ellis SE, Sokolowski DJ, Davis S, Wilson MD, Leek JT, Jaffe
#> AE, Collado-Torres L (2019). "recount-brain: a curated repository of
#> human brain RNA-seq datasets metadata." _bioRxiv_. doi: 10.1101/618025
#> (URL: https://doi.org/10.1101/618025), <URL:
#> https://doi.org/10.1101/618025>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {recount-brain: a curated repository of human brain RNA-seq datasets metadata},
#>     author = {Ashkaun Razmara and Shannon E. Ellis and Dustin J. Sokolowski and Sean Davis and Michael D. Wilson and Jeffrey T. Leek and Andrew E. Jaffe and Leonardo Collado-Torres},
#>     year = {2019},
#>     journal = {bioRxiv},
#>     doi = {10.1101/618025},
#>     url = {https://doi.org/10.1101/618025},
#>   }
#> 
#> Imada E, Sanchez DF, Collado-Torres L, Wilks C, Matam T, Dinalankara W,
#> Stupnikov A, Lobo-Pereira F, Yip C, Yasuzawa K, Kondo N, Itoh M, Suzuki
#> H, Kasukawa T, Hon CC, de Hoon MJ, Shin JW, Carninci P, Jaffe AE, Leek
#> JT, Favorov A, Franco GR, Langmead B, Marchionni L (2020). "Recounting
#> the FANTOM CAGE–Associated Transcriptome." _Genome Research_. doi:
#> 10.1101/gr.254656.119 (URL: https://doi.org/10.1101/gr.254656.119),
#> <URL: https://doi.org/10.1101/gr.254656.119>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Recounting the FANTOM CAGE–Associated Transcriptome},
#>     author = {Eddie-Luidy Imada and Diego Fernando Sanchez and Leonardo Collado-Torres and Christopher Wilks and Tejasvi Matam and Wikum Dinalankara and Aleksey Stupnikov and Francisco Lobo-Pereira and Chi-Wai Yip and Kayoko Yasuzawa and Naoto Kondo and Masayoshi Itoh and Harukazu Suzuki and Takeya Kasukawa and Chung Chau Hon and Michiel JL {de Hoon} and Jay W Shin and Piero Carninci and Andrew E. Jaffe and Jeffrey T. Leek and Alexander Favorov and Glória R Franco and Ben Langmead and Luigi Marchionni},
#>     year = {2020},
#>     journal = {Genome Research},
#>     doi = {10.1101/gr.254656.119},
#>     url = {https://doi.org/10.1101/gr.254656.119},
#>   }
```

Please note that the `recount` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the recount project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Development tools

  - Continuous code testing is possible thanks to [GitHub
    actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
    through *[usethis](https://CRAN.R-project.org/package=usethis)*,
    *[remotes](https://CRAN.R-project.org/package=remotes)*,
    *[sysreqs](https://github.com/r-hub/sysreqs)* and
    *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)*
    customized to use [Bioconductor’s docker
    containers](https://www.bioconductor.org/help/docker/) and
    *[BiocCheck](https://bioconductor.org/packages/3.11/BiocCheck)*.
  - Code coverage assessment is possible thanks to
    [codecov](https://codecov.io/gh) and
    *[covr](https://CRAN.R-project.org/package=covr)*.
  - The [documentation website](http://leekgroup.github.io/recount) is
    automatically updated thanks to
    *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
  - The code is styled automatically thanks to
    *[styler](https://CRAN.R-project.org/package=styler)*.
  - The documentation is formatted thanks to
    *[devtools](https://CRAN.R-project.org/package=devtools)* and
    *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

## Teams involved

  - [Jeff Leek’s lab at JHBSPH Biostatistics
    Department](http://jtleek.com/),
  - [Ben Langmead’s lab at JHU Computer
    Science](http://www.langmead-lab.org/),
  - [Kasper Daniel Hansen’s lab at JHBSPH Biostatistics
    Department](https://www.hansenlab.org/),
  - [Leonardo Collado-Torres](http://lcolladotor.github.io/) and [Andrew
    E. Jaffe](http://aejaffe.com/) from [LIBD](https://www.libd.org/),
  - [Abhinav Nellore’s lab at OHSU](http://nellore.bio/),
  - Data hosted by [SciServer at JHU](https://www.sciserver.org/).

|                                                                                                                                                                               |                                                                                                      |                                                                                                                                                                         |                                                                                                                                                   |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| <a href="http://www.langmead-lab.org/"><img src="http://www.langmead-lab.org/wp-content/uploads/2014/01/Screen-Shot-2014-02-02-at-5.20.13-PM-1024x199.png" width="250px"></a> | <a href="https://www.libd.org/"><img src="http://aejaffe.com/media/LIBD_logo.jpg" width="250px"></a> | <a href="http://nellore.bio/"><img src="https://seekvectorlogo.net/wp-content/uploads/2018/08/oregon-health-science-university-ohsu-vector-logo.png" width="250px"></a> | <a href="https://www.sciserver.org/"><img src="https://skyserver.sdss.org/dr14/en/images/sciserver_logo_inverted_vertical.png" width="250px"></a> |

<!-- Global site tag (gtag.js) - Google Analytics -->

<script async src="https://www.googletagmanager.com/gtag/js?id=UA-78422749-1"></script>

<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'UA-78422749-1');
</script>
