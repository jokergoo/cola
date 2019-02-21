\name{get_stats-ConsensusPartitionList-method}
\alias{get_stats,ConsensusPartitionList-method}
\title{
Get statistics for consensus partitions from all methods
}
\description{
Get statistics for consensus partitions from all methods
}
\usage{
\S4method{get_stats}{ConsensusPartitionList}(object, k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{k}{number of partitions. The value can only be a single value.}

}
\value{
A matrix of partition statistics for a selected k. Rows in the 
matrix correspond to combinations of top-value methods and partition methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
get_stats(cola_rl, k = 2)
}
