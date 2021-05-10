\name{get_matrix-ConsensusPartition-method}
\alias{get_matrix,ConsensusPartition-method}
\title{
Get the original matrix
}
\description{
Get the original matrix
}
\usage{
\S4method{get_matrix}{ConsensusPartition}(object, full = FALSE, include_all_rows = FALSE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{full}{Whether to extract the complete original matrix.}
  \item{include_all_rows}{Internally used.}

}
\value{
A numeric matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
obj = golub_cola["ATC", "skmeans"]
get_matrix(obj)
}
