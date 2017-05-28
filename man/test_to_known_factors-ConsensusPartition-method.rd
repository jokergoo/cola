\name{test_to_known_factors-ConsensusPartition-method}
\alias{test_to_known_factors,ConsensusPartition-method}
\title{
Test correspondance between predicted and known classes
}
\description{
Test correspondance between predicted and known classes
}
\usage{
\S4method{test_to_known_factors}{ConsensusPartition}(object, k, known = object@known_anno)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions}
  \item{known}{a vector or a data frame with known factors}

}
\value{
A matrix of p-values
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
