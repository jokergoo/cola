\name{collect_classes-ConsensusPartitionList-method}
\alias{collect_classes,ConsensusPartitionList-method}
\title{
Collect classes from ConsensusPartitionList object
}
\description{
Collect classes from ConsensusPartitionList object
}
\usage{
\S4method{collect_classes}{ConsensusPartitionList}(object, k,
    top_method = object@top_method, partition_method = object@partition_method, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object returned by \code{\link{run_all_consensus_partition_methods}}.}
  \item{k}{number of partitions}
  \item{top_method}{a vector of top methods}
  \item{partition_method}{a vector of partition methods}
  \item{...}{other arguments.}

}
\value{
No value is returned
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
