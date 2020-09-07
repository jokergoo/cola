\name{compare_signatures-ConsensusPartition-method}
\alias{compare_signatures,ConsensusPartition-method}
\title{
Compare Signatures from Different k
}
\description{
Compare Signatures from Different k
}
\usage{
\S4method{compare_signatures}{ConsensusPartition}(object, k = object@k, verbose = interactive(), ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object. }
  \item{k}{Number of subgroups. Value should be a vector.}
  \item{verbose}{Whether to print message.}
  \item{...}{Other arguments passed to \code{\link{get_signatures,ConsensusPartition-method}}.}

}
\details{
It plots an Euler diagram showing the overlap of signatures from different k.
}
\examples{
\donttest{
data(golub_cola)
res = golub_cola["ATC", "skmeans"]
compare_signatures(res)
}
}
