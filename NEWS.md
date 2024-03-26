# recount 1.29.1

BUG FIXES

* Merged a pull request by @hpages which addresses some changes
to GenomicFeatures. See https://github.com/leekgroup/recount/pull/24
for details.

# recount 1.19.2

BUG FIXES

* Fix a bug in `geo_info()` for reading files on Windows where a trailing
`\r` was added to all variables.
* Avoid the `implicit list embedding of S4 objects is deprecated` warning
that was noted at 
https://github.com/leekgroup/recount/runs/3286046827?check_suite_focus=true#step:20:1417.

# recount 1.13.2

SIGNIFICANT USER-VISIBLE CHANGES

* Documentation website is now available at
http://leekgroup.github.io/recount/. It gets updated with every
commit on the master branch (bioc-devel) using GitHub Actions and pkgdown.

# recount 1.13.1

* Mention in `reproduce_ranges()` the link to
https://support.bioconductor.org/p/126148/#126173 which shows how to update
the gene symbols in the RSE objects in recount.

# recount 1.11.14

* Added a `NEWS.md` file to track changes to the package.

# recount 1.11.13


NEW FEATURES

* Added the function `getTPM()` as discussed in
https://support.bioconductor.org/p/124265
and based on Sonali Arora et al
https://www.biorxiv.org/content/10.1101/445601v2.


# recount 1.11.12


BUG FIXES

* Now `geo_characteristics()` can deal with the scenario reported at
https://support.bioconductor.org/p/116480/ by @Jacques.van-Helden.


# recount 1.11.7


SIGNIFICANT USER-VISIBLE CHANGES

* Renamed `.load_install()` as `.load_check()` as this function now only checks
that the package(s) was installed and returns an error if missing. The
error shows the user how to install the package(s) they are missing
instead of installing them automatically. This complies with Marcel
Ramos' request at https://github.com/leekgroup/recount/issues/14.

# recount 1.11.4


NEW FEATURES

* Added the function `download_retry()` based on
http://bioconductor.org/developers/how-to/web-query/ such that
`download_file()` and other recount functions will re-try to download a
file 3 times before giving up. This should help reduce the number of
occasional failed Bioconductor nightly checks.


# recount 1.9.7


SIGNIFICANT USER-VISIBLE CHANGES

* Cleaned up the documentation of `add_metadata()`
and changed the default `source` from
`recount_brain_v1` to `recount_brain_v2`.

# recount 1.9.6


NEW FEATURES

* `citation('recount')[5]` now lists the `recount-brain`
bioRxiv pre-print citation information.


# recount 1.9.5


NEW FEATURES

* We have now released the FANTOM-CAT/recount2 RSE files
which you can access now through
`download_study(type = 'rse-fc')`. See
Imada EL, Sanchez DF, et al, bioRxiv, 2019
https://www.biorxiv.org/content/10.1101/659490v1
for more information.


# recount 1.9.4


BUG FIXES

* I made the example code for `geo_characteristics()` more robust
since currently `rentrez` can occasionally fails.


# recount 1.9.3


NEW FEATURES

* Add ORCID's following changes at
http://bioconductor.org/developers/package-guidelines/#description

# recount 1.9.2


BUG FIXES

* Added the argument `async` to `snaptron_query()` which
can be set to `FALSE` to address the issue reported
at https://github.com/ChristopherWilks/snaptron/issues/11


# recount 1.7.7


BUG FIXES

* Updated `reproduce_ranges()` to match the `URL` change in
Gencode from
ftp://ftp.sanger.ac.uk to ftp://ftp.ebi.ac.uk


# recount 1.7.5


SIGNIFICANT USER-VISIBLE CHANGES

* `add_metadata()` can now download the `recount_brain_v2` data.


# recount 1.7.4


BUG FIXES

* Fix a `NOTE` about `RefManageR`.


# recount 1.7.3


SIGNIFICANT USER-VISIBLE CHANGES

* Use `BiocManager`

# recount 1.7.2


BUG FIXES

* Fix a unit test.


# recount 1.7.1


SIGNIFICANT USER-VISIBLE CHANGES

* `rse_tx` URLs now point to v2 to reflect recent changes by Fu et al.


# recount 1.5.11


BUG FIXES

* Change some examples to dontrun and improve the code that cleans up after
the tests. This should reduce the size of files left in tmp although
they didn't seem too big to begin with.


# recount 1.5.9


SIGNIFICANT USER-VISIBLE CHANGES

* The functions `add_metadata()` and `add_predictions()` now return the
sample metadata or predictions when the `rse` argument is missing.

# recount 1.5.6


NEW FEATURES

* Added the function `add_metadata()` which can be used to append curated
metadata to a recount rse object. Currently, `add_metadata()` only
supports the `recount_brain_v1` data available at
http://lieberinstitute.github.io/recount-brain/
and to be further described in Razmara et al, in prep, 2018.

# recount 1.5.5


BUG FIXES

* Fix doc link in `geo_characteristics()` which affected the Windows build
machines.

# recount 1.5.4


BUG FIXES

* Fix a unit test for `download_study()`, add another test for the versions,
and fix a `NOTE` in `R CMD check`.

# recount 1.5.3


NEW FEATURES

* `download_study()` can now download the transcript counts (`rse_tx.RData`)
files. The transcript estimation is described in Fu et al, 2018.

SIGNIFICANT USER-VISIBLE CHANGES

* `download_study()` now has a version parameter (defaults to 2). This
argument controls which version of the files to download based on the
change on how exons were defined. Version 1 are reduced exons while
version 2 are disjoint exons as described in further detail in the
documentation tab of the recount website
https://jhubiostatistics.shinyapps.io/recount/.
* `recount_url` and the example `rse_gene_SRP009615` have been updated to match
the changes in version 2.


# recount 1.3.13


BUG FIXES

* Changed `reproduce_ranges()` since disjoint exons are more useful than
reduced exons for downstream analyses.


# recount 1.3.12


NEW FEATURES

* Added the function `read_counts()`.


# recount 1.3.9


SIGNIFICANT USER-VISIBLE CHANGES

* Added citations for
https://www.biorxiv.org/content/early/2017/06/03/145656 and
https://f1000research.com/articles/6-1558/v1 as well as mentions to
them in the vignette.


# recount 1.3.7


SIGNIFICANT USER-VISIBLE CHANGES

* `add_predictions()` was bumped to version 0.0.05


# recount 1.3.5


SIGNIFICANT USER-VISIBLE CHANGES

* Vignette now uses the new `BiocStyle::html_document` that was recently
released.


# recount 1.3.2


NEW FEATURES

* `coverage_matrix()` now has two new arguments: `scale` and `round`. Use
`scale = FALSE` to get raw coverage counts, which you can then scale with
`scale_counts()`. `scale` is set to `TRUE` by default, so the counts are
scaled to a library size of 40 million reads. `round` is set to `FALSE` by
default, but can be set to `TRUE` if you want to get integer counts, just
as in the default of `scale_counts()`.

# recount 1.3.1


SIGNIFICANT USER-VISIBLE CHANGES

* Changed the default version argument of `add_predictions()` to `latest`.
Internally, that's still 0.0.03.


# recount 1.1.27


NEW FEATURES

* Added the `add_predictions()` function which appends the predicted
phenotypes to a RSE object downloaded with recount. The phenotypes
were predicted by Shannon Ellis et al, 2017 (citation coming up soon!).

# recount 1.1.26


SIGNIFICANT USER-VISIBLE CHANGES

* Changed the citation now that the recount2 paper has been published at
http://www.nature.com/nbt/journal/v35/n4/full/nbt.3838.html.


# recount 1.1.25


NEW FEATURES

* Added the function `getRPKM()` which can be used with
`RangedSummarizedExperiment` objects from `recount` and from other sources.

# recount 1.1.24


SIGNIFICANT USER-VISIBLE CHANGES

* `recount_url` now includes the URLs for the GTEx bigWig files.

# recount 1.1.19


SIGNIFICANT USER-VISIBLE CHANGES

* `coverage_matrix()` now returns a RangedSummarizedExperiment object. This
matches the behavior of `recount.bwtool::coverage_matrix_bwtool()` and
is more consistent with the use of RSE objects in recount.

# recount 1.1.18


BUG FIXES

* `coverage_matrix()`'s helper function `.read_pheno()` was failing for some
projects.

# recount 1.1.16


BUG FIXES

* Fixed a bug in the counts in `coverage_matrix()`. They were being
incorrectly multiplied by 100.


# recount 1.1.14


SIGNIFICANT USER-VISIBLE CHANGES

* Completed the change to Gencode v25 annotation for exon and gene counts.


# recount 1.1.13


SIGNIFICANT USER-VISIBLE CHANGES

* We dropped `TxDb.Hsapiens.UCSC.hg38.knownGene` completely from `recount`
and will be using Gencode v25 instead.

# recount 1.1.12


BUG FIXES

* Updated `snaptron_query()` to comply with recent changes in Snaptron.


# recount 1.1.8


SIGNIFICANT USER-VISIBLE CHANGES

* Updated the package so you can now access TCGA data. Now there's over
8 terabytes of data available in the `recount` project!

# recount 1.1.6


SIGNIFICANT USER-VISIBLE CHANGES

* `snaptron_query()` can now access GTEx and TCGA data.

# recount 1.1.5


SIGNIFICANT USER-VISIBLE CHANGES

* Snaptron changed from stingray.cs.jhu.edu:8090 to snaptron.cs.jhu.edu so
`snaptron_query()` has been changed accordingly.


# recount 1.1.2


SIGNIFICANT USER-VISIBLE CHANGES

* The function `reproduce_ranges()` now has the `db` argument. By default
it's set to `TxDb.Hsapiens.UCSC.hg38.knownGene` to reproduce the actual
information used in `recount`. But it can also be used with
`EnsDb.Hsapiens.v79` to use the ENSEMBL annotation. Then with
`coverage_matrix()` you can get the counts for either an updated
`TxDb.Hsapiens.UCSC.hg38.knownGene` or for `EnsDb.Hsapiens.v79` at the
exon and/or gene levels as shown in the vignette.

# recount 1.1.1


SIGNIFICANT USER-VISIBLE CHANGES

* The vignette now describes how to download all the data, how to check
exon-exon junctions by class, and how to use `SciServer` compute
to access all the `recount` data (over 6 TB) via http://www.sciserver.org/


# recount 0.99.30


NEW FEATURES

* Added the function `snaptron_query()` which queries Intropolis via Snaptron
to find if an exon-exon junction is present in the data.


# recount 0.99.29


BUF FIXES

* Fixed an bug in the vignette. Thanks to Michael Love for noticing it!


# recount 0.99.0


NEW FEATURES

* Created the package skeleton for `recount`
* Added the function `reproduce_ranges()` for re-creating the gene or exon
level information used in the `recount` project.
* Added the function `abstract_search()` for identifying SRA projects of
interest by searching the abstracts.
* Added the function `browse_study()` for opening a browser tab for further
exploring a project.
* Added the function `download_study()` for downloading the data from the
`recount` project.
* Added the function `scale_counts()` for properly scaling the counts before
performing a differential expression analysis with the
`RangedSummarizedExperiment` objects hosted in the `recount` project.
* Added the function `expressed_regions()` for defining the expressed regions
in a chromosome for a given SRA study.
* Added the function `coverage_matrix()` for computing the coverage matrix
based on the regions of interest for a given SRA study.
* Added the function `geo_info()` for obtaining sample information from GEO.
* Added the function `find_geo()` for finding the GEO accession id given a
SRA run accession (`id`). This function will be useful for SRA projects
that did not have GEO entries at the time `recount`'s data was created.
* Added the function `geo_characteristics()` for building a `data.frame()` from
`geo_info()`'s results for the characteristics.
* Added the function `all_metadata()` which downloads all the phenotype data
for all projects. This function can be useful for identifying projects
and/or samples of interests.
