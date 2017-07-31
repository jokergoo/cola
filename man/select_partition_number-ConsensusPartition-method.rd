\name{select_partition_number-ConsensusPartition-method}
\alias{select_partition_number,ConsensusPartition-method}
\alias{select_partition_number}
\title{
Several plots for determining the optimized number of partitions
}
\description{
Several plots for determining the optimized number of partitions
}
\usage{
\S4method{select_partition_number}{ConsensusPartition}(object)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}

}
\details{
There are six plots made:

\itemize{
  \item cdf of the consensus matrix under each k
  \item the cophenetic correlation coefficient
  \item PAC score
  \item mean sihouette score
  \item the sum of intra-partition distance
  \item area increase of the area under the cdf of consensus matrix with increasing k
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
