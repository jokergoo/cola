\name{get_best_k-ConsensusPartition-method}
\alias{get_best_k,ConsensusPartition-method}
\title{
Get the best number of partitions
}
\description{
Get the best number of partitions
}
\usage{
\S4method{get_best_k}{ConsensusPartition}(object)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}

}
\details{
It looks for the best k with highest cophenetic correlation coefficient
or lowest PAC score or highest mean silhouette value.
}
\value{
The best k
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
