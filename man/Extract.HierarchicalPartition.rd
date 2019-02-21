\name{[.HierarchicalPartition}
\alias{[.HierarchicalPartition}
\alias{Extract.HierarchicalPartition}
\title{
Subset the HierarchicalPartition object
}
\description{
Subset the HierarchicalPartition object
}
\usage{
\method{[}{HierarchicalPartition}(x, i)
}
\arguments{

  \item{x}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{i}{index. The value should be numeric or a node ID.}

}
\details{
On each node, there is a \code{\link{ConsensusPartition-class}} object.

Note you cannot get a sub-hierarchy of the partition.
}
\value{
A \code{\link{ConsensusPartition-class}} object.
}
\examples{
data(cola_rh)
cola_rh["01"]
cola_rh["2"]
}
