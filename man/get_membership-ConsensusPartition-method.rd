\name{get_membership-ConsensusPartition-method}
\alias{get_membership,ConsensusPartition-method}
\title{
Get membership matrix
}
\description{
Get membership matrix
}
\usage{
\S4method{get_membership}{ConsensusPartition}(object, k, each = FALSE)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}
  \item{each}{whether return the percentage membership matrix which is summarized from all repetitive partitionings or the individual membership in every random partitioning.}

}
\details{
If \code{each == TRUE}, the value in the membership matrix is the probability
to be in one subgroup, while if \code{each == FALSE}, the membership matrix contains the 
class labels for the partitions in all repetitive partitionings with randomly sampling subset
of rows in the matrix.

The percent membership matrix is calculated by \code{\link[clue]{cl_consensus}}.
}
\value{
A membership matrix where rows correspond to the samples in the original matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
get_membership(obj, k = 2)
get_membership(obj, k = 2, each = TRUE)
}
