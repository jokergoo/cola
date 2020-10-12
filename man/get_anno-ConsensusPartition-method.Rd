\name{get_anno-ConsensusPartition-method}
\alias{get_anno,ConsensusPartition-method}
\title{
Get annotations
}
\description{
Get annotations
}
\usage{
\S4method{get_anno}{ConsensusPartition}(object, full = FALSE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{full}{Whether to extract the complete original matrix.}

}
\value{
A data frame if \code{anno} was specified in \code{\link{run_all_consensus_partition_methods}} or \code{\link{consensus_partition}}, or else \code{NULL}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
