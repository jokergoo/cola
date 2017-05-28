\name{get_membership-ConsensusPartition-method}
\alias{get_membership,ConsensusPartition-method}
\alias{get_membership}
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

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions}
  \item{each}{whether return the percentage membership matrix or the membership in every random round}

}
\details{
If \code{each == TRUE}, the value in the membership matrix is the probability
to be in one subgroup, while if \code{each == FALSE}, the membership matrix contains the 
class labels for the partitions in all randomizations.
}
\value{
A membership matrix where rows correspond to the columns in the original matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
