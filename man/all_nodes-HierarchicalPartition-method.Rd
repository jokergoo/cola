\name{all_nodes-HierarchicalPartition-method}
\alias{all_nodes,HierarchicalPartition-method}
\alias{all_nodes}
\title{
All nodes in the hierarchy
}
\description{
All nodes in the hierarchy
}
\usage{
\S4method{all_nodes}{HierarchicalPartition}(object, filter_node = filter_node_param())
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{Depth in the hierarchy.}
  \item{min_n_signatures}{Minimal number of signatures on the node.}

}
\value{
A vector of node ID.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
all_nodes(golub_cola_rh)
}
