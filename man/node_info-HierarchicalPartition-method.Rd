\name{node_info-HierarchicalPartition-method}
\alias{node_info,HierarchicalPartition-method}
\alias{node_info}
\title{
Information on the nodes
}
\description{
Information on the nodes
}
\usage{
\S4method{node_info}{HierarchicalPartition}(object)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}

}
\details{
It returns the following node-level information:

\describe{
  \item{id}{Node id.}
  \item{n_columns}{Number of columns.}
  \item{n_signatures}{Number of signatures.}
  \item{p_signatures}{Percent of signatures.}
  \item{is_leaf}{Whether the node is a leaf}
}
}
\examples{
# There is no example
NULL

}
