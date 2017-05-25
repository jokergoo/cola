\name{ConsensusPartition-class}
\docType{class}
\alias{ConsensusPartitionList}
\alias{ConsensusPartition-class}
\title{
The ConsensusPartition class
}
\description{
The ConsensusPartition class
}
\section{Methods}{
The \code{\link{ConsensusPartition-class}} has following methods:

\describe{
  \item{\code{\link{consensus_partition}}:}{constructor method;}
  \item{\code{\link{plot_ecdf,ConsensusPartition-method}}:}{plot the ecdf of the consensus matrix;}
  \item{\code{\link{select_partition_number,ConsensusPartition-method}}:}{a list of plots used to pick optimized number of partitions;}
  \item{\code{\link{consensus_heatmap,ConsensusPartition-method}}:}{Heatmap of the consensus matrix;}
  \item{\code{\link{membership_heatmap,ConsensusPartition-method}}:}{Heatmap of the membership in each randomization;}
  \item{\code{\link{get_signatures,ConsensusPartition-method}}:}{Heatmap of signature rows;}
  \item{\code{\link{dimension_reduction,ConsensusPartition-method}}:}{dimension reduction plots;}
  \item{\code{\link{collect_plots,ConsensusPartition-method}}:}{Heatmaps for consensus matrix and membership matrix with different number of partitions;}
  \item{\code{\link{collect_classes,ConsensusPartition-method}}:}{Heatmap of classes with different number of partitions;}
  \item{\code{\link{get_param,ConsensusPartition-method}}:}{get parameters for the consensus clustering;}
  \item{\code{\link{get_consensus,ConsensusPartition-method}}:}{get the consensus matrix;}
  \item{\code{\link{get_membership,ConsensusPartition-method}}:}{get the membership in randomizations;}
  \item{\code{\link{get_stat,ConsensusPartition-method}}:}{get statistics for the consensus clustering;}
  \item{\code{\link{get_class,ConsensusPartition-method}}:}{get the consensus class IDs and other columns.}
}}
\examples{
# There is no example
NULL

}
