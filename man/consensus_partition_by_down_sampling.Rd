\name{consensus_partition_by_down_sampling}
\alias{consensus_partition_by_down_sampling}
\title{
Consensus partitioning only with a subset of columns
}
\description{
Consensus partitioning only with a subset of columns
}
\usage{
consensus_partition_by_down_sampling(data,
    top_value_method = "ATC",
    top_n = NULL,
    partition_method = "skmeans",
    max_k = 6, k = NULL,
    subset = min(round(ncol(data)*0.2), 250), pre_select = TRUE,
    verbose = TRUE, prefix = "", anno = NULL, anno_col = NULL,
    dist_method = c("euclidean", "correlation", "cosine"),
    .env = NULL, .predict = TRUE, mc.cores = 1, cores = mc.cores, ...)
}
\arguments{

  \item{data}{A numeric matrix where subgroups are found by columns.}
  \item{top_value_method}{A single top-value method. Available methods are in \code{\link{all_top_value_methods}}. Use \code{\link{register_top_value_methods}} to add a new top-value method.}
  \item{top_n}{Number of rows with top values. The value can be a vector with length > 1. When n > 5000,  the function only randomly sample 5000 rows from top n rows. If \code{top_n} is a vector, paritition will be applied to every values in \code{top_n} and consensus partition is summarized from all partitions.}
  \item{partition_method}{A single partitioning method. Available methods are in \code{\link{all_partition_methods}}. Use \code{\link{register_partition_methods}} to add a new partition method.}
  \item{max_k}{Maximal number of subgroups to try. The function will try for \code{2:max_k} subgroups}
  \item{k}{Alternatively, you can specify a vector k.}
  \item{subset}{Number of columns to randomly sample, or a vector of selected indices.}
  \item{pre_select}{Whether to pre-select by k-means.}
  \item{verbose}{Whether to print messages.}
  \item{prefix}{Internally used.}
  \item{anno}{Annotation data frame.}
  \item{anno_col}{Annotation colors.}
  \item{dist_method}{Method for predict the class for other columns.}
  \item{.env}{An environment, internally used.}
  \item{.predict}{Internally used.}
  \item{mc.cores}{Number of cores. This argument will be removed in future versions.}
  \item{cores}{Number of cores, or a \code{cluster} object returned by \code{\link[parallel]{makeCluster}}.}
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
