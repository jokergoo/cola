\name{all_leaves-HierarchicalPartition-method}
\alias{all_leaves,HierarchicalPartition-method}
\alias{all_leaves}
\title{
All leaves in the hierarchy
}
\description{
All leaves in the hierarchy
}
\usage{
\S4method{all_leaves}{HierarchicalPartition}(object, depth = max_depth(object))
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth in the hierarchy.}

}
\value{
A vector of node ID.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
all_leaves(cola_rh)
}
