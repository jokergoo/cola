\name{get_matrix-ConsensusPartition-method}
\alias{get_matrix,ConsensusPartition-method}
\title{
Get original matrix
}
\description{
Get original matrix
}
\usage{
\S4method{get_matrix}{ConsensusPartition}(object)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}

}
\value{
A numeric matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
obj = cola_rl["sd", "kmeans"]
get_matrix(obj)
}
