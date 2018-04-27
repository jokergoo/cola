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
\S4method{all_nodes}{HierarchicalPartition}(object, depth = NULL)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth in the hierarchy.}

}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
all_nodes(cola_rh)
}
