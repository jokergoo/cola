\name{consensus_partition_by_down_sampling}
\alias{consensus_partition_by_down_sampling}
\title{
Consensus partitioning only with a subset of columns
}
\description{
Consensus partitioning only with a subset of columns
}
\usage{
consensus_partition_by_down_sampling(data, subset = min(round(ncol(data)*0.2), 250),
    verbose = TRUE, prefix = "", anno = NULL, anno_col = NULL,
    dist_method = c("euclidean", "correlation", "cosine"),
    .env = NULL, ...)
}
\arguments{

  \item{data}{A numeric matrix where subgroups are found by columns.}
  \item{subset}{Number of columns to randomly sample.}
  \item{verbose}{Whether to print messages.}
  \item{prefix}{Internally used.}
  \item{anno}{Annotation data frame.}
  \item{anno_col}{Annotation colors.}
  \item{dist_method}{Method for predict the class for other columns.}
  \item{.env}{An environment, internally used.}
  \item{...}{All pass to \code{\link{consensus_partition}}.}

}
\details{
The function performs consensus partitioning only with a small subset
of columns and the class of other columns are predicted by \code{\link{predict_classes,ConsensusPartition-method}}.
}
\examples{
\dontrun{
data(golub_cola)
m = get_matrix(golub_cola)

set.seed(123)
golub_cola_ds = consensus_partition_by_down_sampling(m, subset = 50,
	anno = get_anno(golub_cola), anno_col = get_anno_col(golub_cola),
	top_value_method = "SD", partition_method = "kmeans")
}
}
