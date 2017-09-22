\name{get_stat-ConsensusPartitionList-method}
\alias{get_stat,ConsensusPartitionList-method}
\title{
Get statistics
}
\description{
Get statistics
}
\usage{
\S4method{get_stat}{ConsensusPartitionList}(object, k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{k}{number of partitions.}

}
\value{
A matrix of partition statistics for a selected k. Rows in the 
matrix correspond to all combinations of top methods and partition methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
