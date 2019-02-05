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
    top_value_method = object@top_value_method,
    partition_method = object@partition_method)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object returned by \code{\link{run_all_consensus_partition_methods}}.}
  \item{k}{number of partitions}
  \item{top_value_method}{a vector of top value methods}
  \item{partition_method}{a vector of partition methods}

}
\details{
There are following panels in the plot:

\itemize{
  \item a heatmap shows partitions predicted from all methods where the top annotation is the consensus partition summarized from partitions from all methods, weighted by mean silhouette scores.
  \item a row barplot annotation showing the mean silhouette scores for different methods.
  \item a heatmap shows the similarities of the partitions of pairwise methods, calculated by \code{\link[clue]{cl_dissimilarity}} with the \code{comembership} method.
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rl)
collect_classes(cola_rl, k = 3)
}
}
