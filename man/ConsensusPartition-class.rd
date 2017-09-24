\name{ConsensusPartition-class}
\docType{class}
\alias{ConsensusPartition}
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
  \item{\code{\link{consensus_partition}}:}{constructor method, run consensus partition with a specified top method and a partition method.}
  \item{\code{\link{select_partition_number,ConsensusPartition-method}}:}{make a list of plots used to select optimized number of partitions.}
  \item{\code{\link{consensus_heatmap,ConsensusPartition-method}}:}{Heatmap of the consensus matrix.}
  \item{\code{\link{membership_heatmap,ConsensusPartition-method}}:}{Heatmap of the membership in each random sampling.}
  \item{\code{\link{get_signatures,ConsensusPartition-method}}:}{get the signature rows and make heatmaps.}
  \item{\code{\link{dimension_reduction,ConsensusPartition-method}}:}{dimension reduction plots.}
  \item{\code{\link{collect_plots,ConsensusPartition-method}}:}{Heatmaps for consensus matrix and membership matrix with different number of partitions.}
  \item{\code{\link{collect_classes,ConsensusPartition-method}}:}{Heatmap of classes with different number of partitions.}
  \item{\code{\link{get_param,ConsensusPartition-method}}:}{get parameters for the consensus clustering.}
  \item{\code{\link{get_consensus,ConsensusPartition-method}}:}{get the consensus matrix.}
  \item{\code{\link{get_membership,ConsensusPartition-method}}:}{get the membership in random samplings.}
  \item{\code{\link{get_stat,ConsensusPartition-method}}:}{get statistics for the consensus clustering.}
  \item{\code{\link{get_class,ConsensusPartition-method}}:}{get the consensus class IDs and other columns.}
  \item{\code{\link{get_best_k,ConsensusPartition-method}}:}{guess the best number of partitions.}
  \item{\code{\link{signature_density,ConsensusPartition-method}}:}{plot the density distribution of signatures in different groups.}
  \item{\code{\link{test_to_known_factors,ConsensusPartition-method}}:}{test correlation between predicted subgrouping and known annotations, if available.}
}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
