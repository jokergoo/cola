\name{is_stable_k-ConsensusPartition-method}
\alias{is_stable_k,ConsensusPartition-method}
\title{
Test whether the current k corresponds to a stable partition
}
\description{
Test whether the current k corresponds to a stable partition
}
\usage{
\S4method{is_stable_k}{ConsensusPartition}(object, k, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{...}{Pass to \code{\link{suggest_best_k,ConsensusPartition-method}}.}

}
\details{
if 1-PAC for the k is larger than 0.9 (10\% ambiguity for
the partition), cola marks it as a stable partition.
}
\value{
Logical scalar.
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
is_stable_k(obj, k = 2)
is_stable_k(obj, k = 3)
}
