\name{suggest_best_k-ConsensusPartitionList-method}
\alias{suggest_best_k,ConsensusPartitionList-method}
\title{
Suggest the best number of partitions
}
\description{
Suggest the best number of partitions
}
\usage{
\S4method{suggest_best_k}{ConsensusPartitionList}(object, rand_index_cutoff = 0.95)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{rand_index_cutoff}{The cutoff for Rand index compared to previous k.}

}
\details{
It basically gives the best k for each combination of top-value method and partition method by calling \code{\link{suggest_best_k,ConsensusPartition-method}}.

1-PAC score higher than 0.95 is treated as very stable partition and higher than 0.9 is treated as stable partition.
}
\value{
A data frame with the best k and other statistics for each combination of methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
suggest_best_k(cola_rl)
}
