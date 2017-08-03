\name{get_single_run-ConsensusPartitionList-method}
\alias{get_single_run,ConsensusPartitionList-method}
\title{
Get result for a single top method and partition method
}
\description{
Get result for a single top method and partition method
}
\usage{
\S4method{get_single_run}{ConsensusPartitionList}(object, top_method = object@top_method[1],
    partition_method = object@partition_method[1])
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{top_method}{a single string which is used in \code{\link{run_all_consensus_partition_methods}}}
  \item{partition_method}{a single string which is used in \code{\link{run_all_consensus_partition_methods}}}

}
\value{
A \code{\link{ConsensusPartition-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
