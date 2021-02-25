\name{[[.HierarchicalPartition}
\alias{ExtractExtract.HierarchicalPartition}
\alias{[[.HierarchicalPartition}
\title{
Subset the HierarchicalPartition object
}
\description{
Subset the HierarchicalPartition object
}
\usage{
\method{[[}{HierarchicalPartition}(x, i)
}
\arguments{

  \item{x}{A \code{\link{HierarchicalPartition-class}} object}
  \item{i}{Index. The value should be numeric or a node ID.}

}
\details{
On each node, there is a \code{\link{ConsensusPartition-class}} object.

Note you cannot get a sub-hierarchy of the partition.
}
\value{
A \code{\link{ConsensusPartition-class}} object.
}
\examples{
# There is no example
NULL

}
