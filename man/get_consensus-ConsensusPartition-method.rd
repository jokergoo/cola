\name{get_consensus-ConsensusPartition-method}
\alias{get_consensus,ConsensusPartition-method}
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

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}

}
\details{
For row i and column j in the consensus matrix, the value of corresponding x_ij
is the probability of sample i and sample j being in a same group from all partitions.
}
\value{
A consensus matrix corresponding to the current k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
get_consensus(obj, k = 2)
}
\alias{get_consensus}
