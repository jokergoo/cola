\name{compare_partitions-ConsensusPartition-method}
\alias{compare_partitions,ConsensusPartition-method}
\alias{compare_partitions}
\title{
Compare two partitionings
}
\description{
Compare two partitionings
}
\usage{
\S4method{compare_partitions}{ConsensusPartition}(object, object2, output_file, k1 = 2, k2 = 2,
    dimension_reduction_method = "UMAP",
    id_mapping = guess_id_mapping(rownames(object), "org.Hs.eg.db", FALSE),
    row_km1 = ifelse(k1 == 2, 2, 1),
    row_km2 = ifelse(k1 ==2 && k2 == 2, 2, 1),
    row_km3 = ifelse(k2 == 2, 2, 1))
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition}} object.}
  \item{object2}{A \code{\link{ConsensusPartition}} object.}
  \item{output_file}{The path of the output HTML file. If it is not specified, the report will be opened in the web browser.}
  \item{k1}{Number of subgroups in \code{object}.}
  \item{k2}{Number of subgroups in \code{object2}.}
  \item{dimension_reduction_method}{Which dimension reduction method to use.}
  \item{id_mapping}{Pass to \code{\link{functional_enrichment,ConsensusPartition-method}}.}
  \item{row_km1}{Number of k-means groups, see Details.}
  \item{row_km2}{Number of k-means groups, see Details.}
  \item{row_km3}{Number of k-means groups, see Details.}

}
\details{
The function produces a HTML report which includes comparisons between two partitioning results.

In the report, there are three heatmaps which visualize A) the signature genes specific in the first partition, B) the signature genes
both in the two partitionings and C) the signatures genes specific in the second partition. Argument \code{row_km1}, \code{row_km2} and 
\code{row_km3} control how many k-means groups should be applied on the three heatmaps.
}
\examples{
\dontrun{
data(golub_cola)
require(hu6800.db)
x = hu6800ENTREZID
mapped_probes = mappedkeys(x)
id_mapping = unlist(as.list(x[mapped_probes]))
compare_partitions(golub_cola["ATC:skmeans"], golub_cola["SD:kmeans"], 
    id_mapping = id_mapping)
}
}
