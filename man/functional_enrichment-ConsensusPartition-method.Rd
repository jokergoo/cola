\name{functional_enrichment-ConsensusPartition-method}
\alias{functional_enrichment,ConsensusPartition-method}
\title{
Perform Gene Ontology Enrichment on Signature Genes
}
\description{
Perform Gene Ontology Enrichment on Signature Genes
}
\usage{
\S4method{functional_enrichment}{ConsensusPartition}(object, gene_fdr_cutoff = cola_opt$fdr_cutoff, k = suggest_best_k(object),
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
  \item{ontology}{"BP": biological processes, "MF": molecular functions, "CC": cellular components. }
  \item{min_set_size}{The minimal size of the GO gene sets.}
  \item{max_set_size}{The maximal size of the GO gene sets.}
  \item{verbose}{Whether to print messages.}
  \item{...}{Other arguments.}

}
\value{
A list of three data frames which correspond to results for three GO catalogues:

\itemize{
  \item \code{BP}: biological processes
  \item \code{MF}: molecular functions
  \item \code{CC}: cellular components
}
}
\examples{
# There is no example
NULL

}
