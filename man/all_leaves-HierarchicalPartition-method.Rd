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
\S4method{all_leaves}{HierarchicalPartition}(object, merge_node = merge_node_param())
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{merge_node}{Parameters to merge sub-dendrograms, see \code{\link{merge_node_param}}.}

}
\value{
A vector of node ID.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_rh)
all_leaves(golub_cola_rh)
}
