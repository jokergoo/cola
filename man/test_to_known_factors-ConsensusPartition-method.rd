\name{test_to_known_factors-ConsensusPartition-method}
\alias{test_to_known_factors,ConsensusPartition-method}
\title{
Test correspondance between predicted and known classes
}
\description{
Test correspondance between predicted and known classes
}
\usage{
\S4method{test_to_known_factors}{ConsensusPartition}(object, k, known = object@known_anno, silhouette_cutoff = 0.5)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions}
  \item{known}{a vector or a data frame with known factors}
  \item{silhouette_cutoff}{cutoff for sihouette scores. Samples with value less than this are omit.}

}
\value{
A data frame with columns:
- the number of samples used to test after filtering by \code{silhouette_cutoff}
- p-values from the tests
- number of partitions
}
\seealso{
\code{\link{test_between_factors}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
