\name{ConsensusPartitionList-class}
\docType{class}
\alias{ConsensusPartition}
\alias{ConsensusPartitionList-class}
\title{
The ConsensusPartitionList class
}
\description{
The ConsensusPartitionList class
}
\details{
The object contains results from all combinations of top methods
and partition methods.
}
\section{Methods}{
The \code{\link{ConsensusPartitionList-class}} provides following methods:

\describe{
  \item{\code{\link{run_all_consensus_partition_methods}}:}{constructor method;}
  \item{\code{\link{get_single_run,ConsensusPartitionList-method}}:}{get a single \code{\link{ConsensusPartition-class}} object from the list;}
  \item{\code{\link{top_rows_overlap,ConsensusPartitionList-method}}:}{plot the overlaps of top rows under different top methods;}
  \item{\code{\link{top_rows_heatmap,ConsensusPartitionList-method}}:}{plot the heatmap of top rows under different top methods;}
  \item{\code{\link{get_class,ConsensusPartitionList-method}}:}{get a consensus class IDs}
  \item{\code{\link{collect_plots,ConsensusPartitionList-method}}:}{collect plots from all combination of top methods and partition methods;}
  \item{\code{\link{collect_classes,ConsensusPartitionList-method}}:}{make a plot which contains predicted classes from all combination of top methods and partition methods.}
}}
\seealso{
The \code{\link{ConsensusPartition-class}}.
}
\examples{
# There is no example
NULL

}
