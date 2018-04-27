\name{get_membership-ConsensusPartitionList-method}
\alias{get_membership,ConsensusPartitionList-method}
\title{
Get membership matrix
}
\description{
Get membership matrix
}
\usage{
\S4method{get_membership}{ConsensusPartitionList}(object, k)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{k}{number of partitions.}

}
\details{
The membership matrix (the probability of each sample to be in a subgroup) is inferred
from the membership matrices of all combinations of methods, weighted by the mean silhouette score of the partitions
for each method. So methods which give instable partitions have lower weights 
when summarizing membership matrix from all methods.
}
\value{
A membership matrix where rows correspond to the samples in the original matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
get_membership(cola_rl, k = 2)
}
