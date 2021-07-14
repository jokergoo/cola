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

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{node_id}{A single ID of a node that is going to be split.}
  \item{subset}{The same as in \code{\link{hierarchical_partition}}.}
  \item{min_samples}{The same as in \code{\link{hierarchical_partition}}.}
  \item{max_k}{max_k The same as in \code{\link{hierarchical_partition}}.}
  \item{cores}{Number of cores.}
  \item{verbose}{Whether to print messages.}
  \item{top_n}{The same as in \code{\link{hierarchical_partition}}.}
  \item{min_n_signatures}{The same as in \code{\link{hierarchical_partition}}.}
  \item{group_diff}{The same as in \code{\link{hierarchical_partition}}.}
  \item{fdr_cutoff}{The same as in \code{\link{hierarchical_partition}}.}

}
\details{
It applies hierarchical consensus partitioning on the specified node.
}
\value{
A \code{\link{HierarchicalPartition-class}} object.
}
\examples{
# There is no example
NULL

}
