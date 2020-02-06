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
  \item{\code{\link{consensus_partition}}:}{constructor method, run consensus partitioning with a specified top-value method and a partitioning method.}
  \item{\code{\link{select_partition_number,ConsensusPartition-method}}:}{make a list of plots for selecting optimized number of subgroups.}
  \item{\code{\link{consensus_heatmap,ConsensusPartition-method}}:}{make heatmap of the consensus matrix.}
  \item{\code{\link{membership_heatmap,ConsensusPartition-method}}:}{make heatmap of the membership for individual partitions.}
  \item{\code{\link{get_signatures,ConsensusPartition-method}}:}{get the signature rows and make heatmap.}
  \item{\code{\link{dimension_reduction,ConsensusPartition-method}}:}{make dimension reduction plots.}
  \item{\code{\link{collect_plots,ConsensusPartition-method}}:}{make heatmaps for consensus matrix and membership matrix with different number of subgroups.}
  \item{\code{\link{collect_classes,ConsensusPartition-method}}:}{make heatmap with subgroups with different numbers.}
  \item{\code{\link{get_param,ConsensusPartition-method}}:}{get parameters for the consensus clustering.}
  \item{\code{\link{get_matrix,ConsensusPartition-method}}:}{get the original matrix.}
  \item{\code{\link{get_consensus,ConsensusPartition-method}}:}{get the consensus matrix.}
  \item{\code{\link{get_membership,ConsensusPartition-method}}:}{get the membership of partitions generated from random samplings.}
  \item{\code{\link{get_stats,ConsensusPartition-method}}:}{get statistics for the consensus partitioning.}
  \item{\code{\link{get_classes,ConsensusPartition-method}}:}{get the consensus subgroup labels and other columns.}
  \item{\code{\link{suggest_best_k,ConsensusPartition-method}}:}{guess the best number of subgroups.}
  \item{\code{\link{test_to_known_factors,ConsensusPartition-method}}:}{test correlation between predicted subgroups and known factors, if available.}
  \item{\code{\link{cola_report,ConsensusPartition-method}}:}{generate a HTML report for the whole analysis.}
  \item{\code{\link{functional_enrichment,ConsensusPartition-method}}:}{perform functional enrichment analysis on significant genes if rows in the matrix can be corresponded to genes.}
}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
