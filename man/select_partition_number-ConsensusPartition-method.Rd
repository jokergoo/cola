\name{select_partition_number-ConsensusPartition-method}
\alias{select_partition_number,ConsensusPartition-method}
\title{
Several plots for determining the optimized number of partitions
}
\description{
Several plots for determining the optimized number of partitions
}
\usage{
\S4method{select_partition_number}{ConsensusPartition}(object, all_stats = FALSE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{all_stats}{Whether to show all statistics that were calculated. Used internally.}

}
\details{
There are following plots made:

\itemize{
  \item ECDF of the consensus matrix under each k, made by \code{\link{plot_ecdf,ConsensusPartition-method}},
  \item \code{\link{PAC}} score,
  \item mean sihouette score,
  \item the \code{\link{concordance}} for each partition to the consensus partition,
  \item area increase of the area under the ECDF of consensus matrix with increasing k,
  \item Rand index for current k compared to k - 1,
  \item Jaccard coefficient for current k compared to k - 1,
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
select_partition_number(cola_rl["sd", "hclust"])
}
\alias{select_partition_number}
