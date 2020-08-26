\name{is_stable_k-ConsensusPartitionList-method}
\alias{is_stable_k,ConsensusPartitionList-method}
\title{
Test whether the current k corresponds to a stable partition
}
\description{
Test whether the current k corresponds to a stable partition
}
\usage{
\S4method{is_stable_k}{ConsensusPartitionList}(object, k, ...)
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
Logical vector
}
\examples{
data(golub_cola)
is_stable_k(golub_cola, k = 3)
}
