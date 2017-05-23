\name{collect_plots-ConsensusPartitionList-method}
\alias{collect_plots,ConsensusPartitionList-method}
\title{
Collect plots from ConsensusPartitionList object
}
\description{
Collect plots from ConsensusPartitionList object
}
\usage{
\S4method{collect_plots}{ConsensusPartitionList}(object, k = 2, fun = consensus_heatmap,
    top_method = object@top_method, partition_method = object@partition_method, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object from \code{\link{run_all_consensus_partition_methods}}.}
  \item{k}{number of partitions.}
  \item{fun}{function used to generate plots. Valid functions are \code{\link{consensus_heatmap,ConsensusPartition-method}}, \code{\link{plot_ecdf,ConsensusPartition-method}}, \code{\link{membership_heatmap,ConsensusPartition-method}} and \code{\link{get_signatures,ConsensusPartition-method}}.}
  \item{top_method}{a vector of top methods.}
  \item{partition_method}{a vector of partition methods.}
  \item{...}{other arguments passed to corresponding \code{fun}.}

}
\examples{
# There is no example
NULL

}
