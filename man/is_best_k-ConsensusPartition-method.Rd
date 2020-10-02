\name{is_best_k-ConsensusPartition-method}
\alias{is_best_k,ConsensusPartition-method}
\title{
Test whether the current k is the best/optional k
}
\description{
Test whether the current k is the best/optional k
}
\usage{
\S4method{is_best_k}{ConsensusPartition}(object, k, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{...}{Pass to \code{\link{suggest_best_k,ConsensusPartition-method}}.}

}
\details{
Optional best k is also assigned as \code{TRUE}.
}
\value{
Logical scalar.
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
is_best_k(obj, k = 2)
is_best_k(obj, k = 3)
}
