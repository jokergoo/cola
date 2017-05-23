\name{get_class-ConsensusPartitionList-method}
\alias{get_class,ConsensusPartitionList-method}
\title{
Get class from the consensus_partition_all_methods object
}
\description{
Get class from the consensus_partition_all_methods object
}
\usage{
\S4method{get_class}{ConsensusPartitionList}(object, k, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{k}{number of partitions}
  \item{...}{other arguments}

}
\details{
The class IDs is re-calculated by merging class IDs from all methods.
}
\value{
A data frame with class IDs and other columns.
}
\examples{
# There is no example
NULL

}
