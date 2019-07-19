\name{get_stats-ConsensusPartition-method}
\alias{get_stats,ConsensusPartition-method}
\title{
Get statistics for the consensus partition
}
\description{
Get statistics for the consensus partition
}
\usage{
\S4method{get_stats}{ConsensusPartition}(object, k = object@k, all_stats = FALSE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of partitions. The value can be a vector.}
  \item{all_stats}{Whether to show all statistics that were calculated. Used internally.}

}
\details{
The statistics are:

\describe{
  \item{PAC}{proportion of ambiguous clustering, calculated by \code{\link{PAC}}.}
  \item{mean_silhouette}{the mean silhouette score. See \url{https://en.wikipedia.org/wiki/Silhouette_(clustering)} .}
  \item{concordance}{the mean probability that each partition fits the consensus partition, calculated by \code{\link{concordance}}.}
  \item{area_increased}{the increased area under ECDF (the empirical cumulative distribution function curve) to the previous k.}
  \item{Rand}{the Rand index which is the percent of pairs of samples that are both in a same cluster or both are not  in a same cluster in the partition of \code{k} and \code{k-1}. See \url{https://en.wikipedia.org/wiki/Rand_index} .}
  \item{Jaccard}{the ratio of pairs of samples are both in a same cluster in the partition of \code{k} and \code{k-1} and the pairs of samples are both in a same cluster in the partition \code{k} or \code{k-1}.}
}
}
\value{
A matrix of partition statistics.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
get_stats(obj)
get_stats(obj, k = 2)
}
