\name{GO_enrichment-ConsensusPartitionList-method}
\alias{GO_enrichment,ConsensusPartitionList-method}
\title{
Perform Gene Ontology Enrichment on Signature Genes
}
\description{
Perform Gene Ontology Enrichment on Signature Genes
}
\usage{
\S4method{GO_enrichment}{ConsensusPartitionList}(object, gene_fdr_cutoff = 0.05,
    id_mapping = guess_id_mapping(rownames(object), org_db, FALSE),
    org_db = "org.Hs.eg.db", ontology = c("BP", "MF", "CC"),
    min_set_size = 10, max_set_size = 1000)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object from \code{\link{run_all_consensus_partition_methods}}.}
  \item{gene_fdr_cutoff}{Cutoff of FDR to define significant signature genes.}
  \item{id_mapping}{If the gene IDs which are row names of the original matrix are not Entrez IDs, a named vector should be provided where the names are the gene IDs in the matrix and values are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.}
  \item{org_db}{Annotation database.}
  \item{ontology}{"BP": biological processes, "MF": molecular functions, "CC": cellular components. }
  \item{min_set_size}{The minimal size of the GO gene sets.}
  \item{max_set_size}{The maximal size of the GO gene sets.}

}
\details{
For each method, the signature genes are extracted based on the best k.

It calls \code{\link{GO_enrichment,ConsensusPartition-method}} on the consensus partitioning results for each method.
}
\value{
A list where each element in the list corresponds to enrichment results from a single method.
}
\examples{
# There is no example
NULL

}
