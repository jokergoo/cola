\name{get_signatures-HierarchicalPartition-method}
\alias{get_signatures,HierarchicalPartition-method}
\title{
Get signatures rows
}
\description{
Get signatures rows
}
\usage{
\S4method{get_signatures}{HierarchicalPartition}(object, depth = NULL)
}
\arguments{

  \item{object}{a \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{minimal depth of the hierarchy}

}
\details{
The function called \code{\link{get_signatures,ConsensusPartition-method}} to find signatures at
each level of the partition hierarchy.
}
\value{
A list of signature names (row names of the original data matrix)
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
