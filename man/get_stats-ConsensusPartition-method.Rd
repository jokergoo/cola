\name{get_stats-ConsensusPartition-method}
\alias{get_stats,ConsensusPartition-method}
\title{
Get statistics
}
\description{
Get statistics
}
\usage{
\S4method{get_stats}{ConsensusPartition}(object, k = object@k, all_stats = FALSE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups. The value can be a vector.}
  \item{all_stats}{Whether to show all statistics that were calculated. Used internally.}

}
\details{
The statistics are:

\describe{
  \item{1-PAC}{1 - proportion of ambiguous clustering, calculated by \code{\link{PAC}}.}
  \item{mean_silhouette}{The mean silhouette score. See \url{https://en.wikipedia.org/wiki/Silhouette_(clustering)} .}
  \item{concordance}{The mean probability that each partition fits the consensus partition, calculated by \code{\link{concordance}}.}
  \item{area_increased}{The increased area under eCDF (the empirical cumulative distribution function) curve to the previous k.}
  \item{Rand}{This is the percent of pairs of samples that are both in a same cluster or both are not  in a same cluster in the partition of \code{k} and \code{k-1}. See \url{https://en.wikipedia.org/wiki/Rand_index} .}
  \item{Jaccard}{The ratio of pairs of samples are both in a same cluster in the partition of \code{k} and \code{k-1} and the pairs of samples are both in a same cluster in the partition \code{k} or \code{k-1}.}
}
}
\value{
A matrix of partition statistics.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
get_stats(obj)
get_stats(obj, k = 2)
}
