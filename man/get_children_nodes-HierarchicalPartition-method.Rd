\name{get_children_nodes-HierarchicalPartition-method}
\alias{get_children_nodes,HierarchicalPartition-method}
\alias{get_children_nodes}
\title{
Test whether a node is a leaf node
}
\description{
Test whether a node is a leaf node
}
\usage{
\S4method{get_children_nodes}{HierarchicalPartition}(object, node, merge_node = merge_node_param())
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{node}{A vector of node IDs.}
  \item{merge_node}{Parameters to merge sub-dendrograms, see \code{\link{merge_node_param}}.}

}
\value{
A vector of children nodes.
}
\examples{
# There is no example
NULL

}
