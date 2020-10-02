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

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{each}{Whether to return the percentage membership matrix which is summarized from all partitions or the individual membership in every single partition run.}

}
\details{
If \code{each == FALSE}, the value in the membership matrix is the probability
to be in one subgroup, while if \code{each == TRUE}, the membership matrix contains the 
subgroup labels for every single partitions which are from randomly sampling from the original matrix.

The percent membership matrix is calculated by \code{\link[clue]{cl_consensus}}.
}
\value{
\itemize{
  \item If \code{each == FALSE}, it returns a membership matrix where rows correspond to the columns from the subgroups.
  \item If \code{each == TRUE}, it returns a membership matrix where rows correspond to the columns from the original matrix.
}
}
\seealso{
\code{\link{get_membership,ConsensusPartitionList-method}} summarizes membership from partitions from all combinations
of top-value methods and partitioning methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
get_membership(obj, k = 2)
get_membership(obj, k = 2, each = TRUE)
}
