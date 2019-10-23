\name{GO_enrichment-HierarchicalPartition-method}
\alias{GO_enrichment,HierarchicalPartition-method}
\title{
Perform Gene Ontology Enrichment on Signature Genes
}
\description{
Perform Gene Ontology Enrichment on Signature Genes
}
\usage{
\S4method{GO_enrichment}{HierarchicalPartition}(object, gene_fdr_cutoff = cola_opt$fdr_cutoff,
    id_mapping = guess_id_mapping(rownames(object), org_db, FALSE),
    org_db = "org.Hs.eg.db",
    min_set_size = 10, max_set_size = 1000)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{gene_fdr_cutoff}{Cutoff of FDR to define significant signature genes.}
  \item{id_mapping}{If the gene IDs which are row names of the original matrix are not Entrez IDs, a named vector should be provided where the names are the gene IDs in the matrix and values are correspoinding Entrez IDs. The value can also be a function that converts gene IDs.}
  \item{org_db}{Annotation database.}
  \item{min_set_size}{The minimal size of the GO gene sets.}
  \item{max_set_size}{The maximal size of the GO gene sets.}

}
\details{
On each node of the partition hierarchy, the signature genes are extracted based on the best k.

It calls \code{\link{GO_enrichment,ConsensusPartition-method}} on the consensus partitioning results on each node.
}
\value{
A list where each element in the list corresponds to enrichment results on a node.
}
\examples{
# There is no example
NULL
}
