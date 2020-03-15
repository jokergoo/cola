\name{compare_signatures-ConsensusPartition-method}
\alias{compare_signatures,ConsensusPartition-method}
\alias{compare_signatures}
\title{
Compare Signatures from Different k
}
\description{
Compare Signatures from Different k
}
\usage{
\S4method{compare_signatures}{ConsensusPartition}(object, k = object@k, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object. }
  \item{k}{Number of subgroups. Value should be a vector.}
  \item{...}{Other arguments passed to \code{\link{get_signatures,ConsensusPartition-method}}.}

}
\details{
It plots an Euler diagram showing the overlap of signatures from different k.
}
\examples{
# There is no example
NULL

}
