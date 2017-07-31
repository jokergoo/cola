\name{test_to_known_factors-ConsensusPartitionList-method}
\alias{test_to_known_factors,ConsensusPartitionList-method}
\title{
Test correspondance between predicted and known classes
}
\description{
Test correspondance between predicted and known classes
}
\usage{
\S4method{test_to_known_factors}{ConsensusPartitionList}(object, k, known = object@list[[1]]@known_anno,
    silhouette_cutoff = 0.5)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{k}{number of partitions}
  \item{known}{a vector or a data frame with known factors}
  \item{silhouette_cutoff}{cutoff for sihouette scores}

}
\details{
The function basically sends each \code{\link{ConsensusPartition-class}} object to
\code{\link{test_to_known_factors,ConsensusPartition-method}} and merges afterwards.
}
\value{
A matrix of p-values where the first column is the number of data columns
to test after filtering by \code{silhouette_cutoff}.
}
\seealso{
\code{\link{test_between_factors}}, \code{\link{test_to_known_factors,ConsensusPartition-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
