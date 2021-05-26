\name{split_node-HierarchicalPartition-method}
\alias{split_node,HierarchicalPartition-method}
\alias{split_node}
\title{
Split node
}
\description{
Split node
}
\usage{
\S4method{split_node}{HierarchicalPartition}(object, node_id,
    subset = object@param$subset,
    min_samples = object@param$min_samples, max_k = object@param$max_k, cores = object@param$cores,
    verbose = TRUE,
    top_n = object@param$top_n, min_n_signatures = object@param$min_n_signatures,
    group_diff = object@param$group_diff, fdr_cutoff = object@param$fdr_cutoff)
}
\arguments{

  \item{object}{object}
  \item{node_id}{node id}
  \item{subset}{subset}
  \item{min_samples}{-min_samples}
  \item{max_k}{max_k}
  \item{cores}{cores}
  \item{verbose}{verobse}
  \item{top_n}{top_n}
  \item{min_n_signatures}{min-sigatures}
  \item{group_diff}{group diff}
  \item{fdr_cutoff}{fdr cutoff}

}
\examples{
# There is no example
NULL

}
