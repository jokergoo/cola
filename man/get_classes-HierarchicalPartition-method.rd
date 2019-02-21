\name{get_classes-HierarchicalPartition-method}
\alias{get_classes,HierarchicalPartition-method}
\title{
Get class IDs from the HierarchicalPartition object
}
\description{
Get class IDs from the HierarchicalPartition object
}
\usage{
\S4method{get_classes}{HierarchicalPartition}(object, depth = max_depth(object))
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{depth of the hierarchy.}

}
\value{
A vector of classes IDs. The class IDs are the node IDs where the subgroup sits in the hierarchy.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
get_classes(cola_rh)
get_classes(cola_rh, depth = 2)
}
