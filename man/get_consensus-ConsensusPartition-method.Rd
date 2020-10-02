\name{get_consensus-ConsensusPartition-method}
\alias{get_consensus,ConsensusPartition-method}
\alias{get_consensus}
\title{
Get consensus matrix
}
\description{
Get consensus matrix
}
\usage{
\S4method{get_consensus}{ConsensusPartition}(object, k)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}

}
\details{
For row i and column j in the consensus matrix, the value of corresponding x_ij
is the probability of sample i and sample j being in the same group from all partitions.
}
\value{
A consensus matrix corresponding to the current k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
get_consensus(obj, k = 2)
}
