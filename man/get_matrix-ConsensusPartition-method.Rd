\name{get_matrix-ConsensusPartition-method}
\alias{get_matrix,ConsensusPartition-method}
\title{
Get the original matrix
}
\description{
Get the original matrix
}
\usage{
\S4method{get_matrix}{ConsensusPartition}(object)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}

}
\value{
A numeric matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["SD", "kmeans"]
get_matrix(obj)
}
