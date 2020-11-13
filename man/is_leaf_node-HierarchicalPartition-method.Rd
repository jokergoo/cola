\name{is_leaf_node-HierarchicalPartition-method}
\alias{is_leaf_node,HierarchicalPartition-method}
\alias{is_leaf_node}
\title{
Test whether a node is a leaf node
}
\description{
Test whether a node is a leaf node
}
\usage{
\S4method{is_leaf_node}{HierarchicalPartition}(object, node, filter_node = filter_node_param())
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{node}{A vector of node IDs.}
  \item{depth}{Depth in the hierarchy.}
  \item{min_n_signatures}{Minimal number of signatures on the node.}

}
\examples{
data(golub_cola_rh)
is_leaf_node(golub_cola_rh, all_leaves(golub_cola_rh))
}
