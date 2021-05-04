\name{functional_enrichment-ConsensusPartition-method}
\alias{functional_enrichment,ConsensusPartition-method}
\title{
Perform functional enrichment on signature genes
}
\description{
Perform functional enrichment on signature genes
}
\usage{
\S4method{functional_enrichment}{ConsensusPartition}(object, gene_fdr_cutoff = cola_opt$fdr_cutoff, k = suggest_best_k(object, help = FALSE),
    row_km = NULL, id_mapping = guess_id_mapping(rownames(object), org_db, verbose),
    org_db = "org.Hs.eg.db", ontology = "BP",
    min_set_size = 10, max_set_size = 1000,
    verbose = TRUE, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object from \code{\link{run_all_consensus_partition_methods}}.}
  \item{gene_fdr_cutoff}{Cutoff of FDR to define significant signature genes.}
  \item{k}{Number of subgroups.}
  \item{row_km}{Number of row clusterings by k-means to separate the matrix that only contains signatures.}
  \item{id_mapping}{If the gene IDs which are row names of the original matrix are not Entrez IDs, a named vector should be provided where the names are the gene IDs in the matrix and values are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.}
  \item{org_db}{Annotation database.}
  \item{ontology}{See corresponding argumnet in \code{\link{functional_enrichment,ANY-method}}.}
  \item{min_set_size}{The minimal size of the gene sets.}
  \item{max_set_size}{The maximal size of the gene sets.}
  \item{verbose}{Whether to print messages.}
  \item{...}{Pass to \code{\link{functional_enrichment,ANY-method}}.}

}
\details{
For how to control the parameters of functional enrichment, see help page of \code{\link{functional_enrichment,ANY-method}}.
}
\seealso{
\url{http://bioconductor.org/packages/devel/bioc/vignettes/cola/inst/doc/functional_enrichment.html}
}
\value{
A list of data frames which correspond to results for the functional ontologies:
}
\examples{
# There is no example
NULL

}
