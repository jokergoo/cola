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
  \item{each}{whether return the percentage membership matrix which is summarized from all partitions or the individual membership in every random partition.}

}
\details{
If \code{each == TRUE}, the value in the membership matrix is the probability
to be in one class, while if \code{each == FALSE}, the membership matrix contains the 
class labels for every single partitions which are from randomly sampling subset
of rows in the matrix.

The percent membership matrix is calculated by \code{\link[clue]{cl_consensus}}.
}
\value{
A membership matrix where rows correspond to the columns in the original matrix.
}
\seealso{
\code{\link{get_membership,ConsensusPartitionList-method}} summarizes membership from partitions from all combinations
of top-value methods and partition methods.
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
