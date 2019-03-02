\name{get_classes-ConsensusPartition-method}
\alias{get_classes,ConsensusPartition-method}
\title{
Get class IDs from the ConsensusPartition object
}
\description{
Get class IDs from the ConsensusPartition object
}
\usage{
\S4method{get_classes}{ConsensusPartition}(object, k = object@k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}

}
\value{
A data frame with class IDs and other columns which are entropy of the percent membership matrix
and the silhouette scores which measure the stability of a sample to stay in its group.

If \code{k} is not specified, it returns a data frame with class IDs from every k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
get_classes(obj, k = 2)
get_classes(obj)
}
