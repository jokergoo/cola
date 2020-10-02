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
\S4method{all_nodes}{HierarchicalPartition}(object, depth = max_depth(object))
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{Depth in the hierarchy.}

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
