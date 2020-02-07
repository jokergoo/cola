\name{get_classes-ConsensusPartition-method}
\alias{get_classes,ConsensusPartition-method}
\title{
Get subgroup labels
}
\description{
Get subgroup labels
}
\usage{
\S4method{get_classes}{ConsensusPartition}(object, k = object@k)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}

}
\value{
A data frame with subgroup labels and other columns which are entropy of the percent membership matrix
and the silhouette scores which measure the stability of a sample to stay in its group.

If \code{k} is not specified, it returns a data frame with subgroup labels from all k.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["SD", "kmeans"]
get_classes(obj, k = 2)
get_classes(obj)
}
