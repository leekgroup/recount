% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_metadata.R
\name{add_metadata}
\alias{add_metadata}
\title{Add additional curated metadata to a recount rse object}
\usage{
add_metadata(rse, source = "recount_brain_v2", is_tcga = FALSE, verbose = TRUE)
}
\arguments{
\item{rse}{A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment-class}
object as downloaded with \link{download_study}. If this argument is
not specified, the function will return the raw metadata table.}

\item{source}{A valid source name. The only supported options at this
moment are \code{recount_brain_v1} and \code{recount_brain_v2}.}

\item{is_tcga}{Set to \code{TRUE} only when \code{rse} is from TCGA.
Otherwise set to \code{FALSE} (default).}

\item{verbose}{If \code{TRUE} it will print a message of where the
predictions file is being downloaded to.}
}
\value{
A \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment-class}
object with the sample metadata columns appended to the \code{colData()}
slot.
}
\description{
This function appends sample metadata information to a
\link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment-class} from the
recount2 project. The sample metadata comes from curated efforts
independent from the original recount2 project. Currently the only
information comes from the recount_brain project described in more detail
at \url{http://lieberinstitute.github.io/recount-brain/}.
}
\details{
For \code{source = "recount_brain_v1"} and
\code{source = "recount_brain_v2"}, the metadata columns are
described at \url{http://lieberinstitute.github.io/recount-brain/}.
Alternatively, you can explore \code{recount_brain_v2} interactively at
\url{https://jhubiostatistics.shinyapps.io/recount-brain/}.

If you use the recount_brain data please cite the Razmara et al.
bioRxiv, 2019 \url{https://www.biorxiv.org/content/10.1101/618025v1}.
A bib file is available via citation('recount').
}
\examples{

## Add the sample metadata to an example rse_gene object
rse_gene <- add_metadata(rse_gene_SRP009615, "recount_brain_v2")

## Explore the metadata
colData(rse_gene)

## For a list of studies present in recount_brain check
## http://lieberinstitute.github.io/recount-brain/.
## recount_brain_v2 includes GTEx and TCGA brain samples in addition to the
## recount_brain_v1 data, plus ontology information.


## Obtain all the recount_brain_v2 metadata if you want to
## explore the metadata manually
recount_brain_v2 <- add_metadata(source = "recount_brain_v2")
}
\references{
Razmara et al, bioRxiv, 2019.
\url{https://www.biorxiv.org/content/10.1101/618025v1}
}
\author{
Leonardo Collado-Torres
}
