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

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions. The value can be a vector.}

}
\details{
The statistics are:

\describe{
  \item{cophcor}{cophenetic correlation coefficient.}
  \item{PAC}{proportion of ambiguous clustering, calculated by \code{\link{PAC}}.}
  \item{mean_silhouette}{the mean silhouette score.}
}
}
\value{
A matrix of partition statistics for all k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
rl = readRDS(system.file("extdata/example.rds", package = "cola"))
obj = rl["sd", "kmeans"]
get_stat(obj)
get_stat(obj, k = 2)
}
