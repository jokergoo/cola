\name{guess_best_k-ConsensusPartition-method}
\alias{guess_best_k,ConsensusPartition-method}
\title{
Guess the best number of partitions
}
\description{
Guess the best number of partitions
}
\usage{
\S4method{guess_best_k}{ConsensusPartition}(object, rand_index_cutoff = 0.9)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{rand_index_cutoff}{the Rand index compared to previous k is larger than this, it is filtered out.}

}
\details{
The best k is voted from 1) which k has the maximum cophcor value, 2) which k has the minimal PAC value,
3) which k has the maximum mean silhouette value and 4) which k has the maximum concordance value.
}
\value{
The best k
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
guess_best_k(obj)
}
