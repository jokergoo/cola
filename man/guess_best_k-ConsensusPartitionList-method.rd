\name{guess_best_k-ConsensusPartitionList-method}
\alias{guess_best_k,ConsensusPartitionList-method}
\title{
Get the best number of partitions
}
\description{
Get the best number of partitions
}
\usage{
\S4method{guess_best_k}{ConsensusPartitionList}(object, rand_index_cutoff = 0.9)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{rand_index_cutoff}{the Rand index compared to previous k is larger than this, it is filtered out.}

}
\details{
It basically gives best k for each combination of top value method and partition method.
}
\value{
A data frame with best k for each combination of methods
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
guess_best_k(cola_rl)
}
