\name{get_classes-HierarchicalPartition-method}
\alias{get_classes,HierarchicalPartition-method}
\title{
Get class IDs from the HierarchicalPartition object
}
\description{
Get class IDs from the HierarchicalPartition object
}
\usage{
\S4method{get_classes}{HierarchicalPartition}(object, merge_node = merge_node_param())
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{merge_node}{Parameters to merge sub-dendrograms, see \code{\link{merge_node_param}}.}

}
\value{
A data frame of classes IDs. The class IDs are the node IDs where the subgroup sits in the hierarchy.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
get_classes(golub_cola_rh)
}
