\name{get_stat-ConsensusPartition-method}
\alias{get_stat,ConsensusPartition-method}
\title{
Get statistics
}
\description{
Get statistics
}
\usage{
\S4method{get_stat}{ConsensusPartition}(object, k = object@k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions. The value can be a vector.}

}
\details{
The statistics are:

\describe{
  \item{cophcor}{cophenetic correlation coefficient. It measures if hierarchical clustering is applied on the consensus matrix, how good it correlates to the consensus matrix itself.}
  \item{PAC}{proportion of ambiguous clustering, calculated by \code{\link{PAC}}.}
  \item{mean_silhouette}{the mean silhouette score.}
  \item{concordance}{the mean probability that each partition fits the consensus partition, calculated by \code{\link{concordance}}.}
  \item{area_increased}{the increased area under ecdf to the previous k}
  \item{Rand}{the Rand index which is the percent of pairs of samples that are both in a same cluster or both are not  in a same cluster in the partition of \code{k} and \code{k-1}.}
  \item{Jaccard}{the ratio of pairs of samples are both in a same cluster in the partition of \code{k} and \code{k-1} and the pairs of samples are both in a same cluster in the partition \code{k} or \code{k-1}.}
}
}
\value{
A matrix of partition statistics for all k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
get_stat(obj)
get_stat(obj, k = 2)
}
