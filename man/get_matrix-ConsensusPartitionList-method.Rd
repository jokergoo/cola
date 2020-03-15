\name{get_matrix-ConsensusPartitionList-method}
\alias{get_matrix,ConsensusPartitionList-method}
\title{
Get the original matrix
}
\description{
Get the original matrix
}
\usage{
\S4method{get_matrix}{ConsensusPartitionList}(object)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}

}
\value{
A numeric matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
get_matrix(cola_rl)
}
