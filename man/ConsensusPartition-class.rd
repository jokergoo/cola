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
  \item{\code{\link{consensus_partition}}:}{constructor method, run consensus partition with a specified top-value method and a partition method.}
  \item{\code{\link{select_partition_number,ConsensusPartition-method}}:}{make a list of plots to select optimized number of partitions.}
  \item{\code{\link{consensus_heatmap,ConsensusPartition-method}}:}{make heatmap of the consensus matrix.}
  \item{\code{\link{membership_heatmap,ConsensusPartition-method}}:}{make heatmap of the membership in every random sampling.}
  \item{\code{\link{get_signatures,ConsensusPartition-method}}:}{get the signature rows and make heatmap.}
  \item{\code{\link{dimension_reduction,ConsensusPartition-method}}:}{make dimension reduction plots.}
  \item{\code{\link{collect_plots,ConsensusPartition-method}}:}{make heatmaps for consensus matrix and membership matrix with different number of partitions.}
  \item{\code{\link{collect_classes,ConsensusPartition-method}}:}{make heatmap of classes with different numbers of partitions.}
  \item{\code{\link{get_param,ConsensusPartition-method}}:}{get parameters for the consensus clustering.}
  \item{\code{\link{get_matrix,ConsensusPartition-method}}:}{get the original matrix.}
  \item{\code{\link{get_consensus,ConsensusPartition-method}}:}{get the consensus matrix.}
  \item{\code{\link{get_membership,ConsensusPartition-method}}:}{get the membership in random samplings.}
  \item{\code{\link{get_stats,ConsensusPartition-method}}:}{get metrics for the consensus clustering.}
  \item{\code{\link{get_classes,ConsensusPartition-method}}:}{get the consensus class IDs and other columns.}
  \item{\code{\link{suggest_best_k,ConsensusPartition-method}}:}{guess the best number of partitions.}
  \item{\code{\link{test_to_known_factors,ConsensusPartition-method}}:}{test correlation between predicted classes and known factors, if available.}
}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
