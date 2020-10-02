\name{is_best_k-ConsensusPartitionList-method}
\alias{is_best_k,ConsensusPartitionList-method}
\title{
Test whether the current k is the best/optional k
}
\description{
Test whether the current k is the best/optional k
}
\usage{
\S4method{is_best_k}{ConsensusPartitionList}(object, k, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{k}{Number of subgroups.}
  \item{...}{Pass to \code{\link{suggest_best_k,ConsensusPartitionList-method}}.}

}
\details{
It tests on the partitions for every method.
}
\value{
Logical vector.
}
\examples{
data(golub_cola)
is_best_k(golub_cola, k = 3)
}
