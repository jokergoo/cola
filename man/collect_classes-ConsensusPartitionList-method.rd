\name{collect_classes-ConsensusPartitionList-method}
\alias{collect_classes,ConsensusPartitionList-method}
\title{
Collect classes from ConsensusPartitionList object
}
\description{
Collect classes from ConsensusPartitionList object
}
\usage{
\S4method{collect_classes}{ConsensusPartitionList}(object, k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object returned by \code{\link{run_all_consensus_partition_methods}}.}
  \item{k}{number of partitions.}

}
\details{
There are following panels in the plot:

\itemize{
  \item a heatmap showing partitions predicted from all methods where the top annotation is the consensus partition summarized from partitions from all methods, weighted by mean silhouette scores.
  \item a row barplot annotation showing the mean silhouette scores for different methods.
}

The row clustering is applied on the dissimilarity matrix calculated by \code{\link[clue]{cl_dissimilarity}} with the \code{comembership} method.

The brightness of the color corresponds to the silhouette scores for the consensus partition in each method.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
collect_classes(cola_rl, k = 3)
}
