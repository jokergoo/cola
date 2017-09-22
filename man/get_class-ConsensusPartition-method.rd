\name{get_class-ConsensusPartition-method}
\alias{get_class,ConsensusPartition-method}
\title{
Get class from the ConsensusPartition object
}
\description{
Get class from the ConsensusPartition object
}
\usage{
\S4method{get_class}{ConsensusPartition}(object, k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions}

}
\value{
A data frame with class IDs and other columns which are entropy of the membership matrix
and the silhouette scores which measure the stability of a sample to stay in its subgroup.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
