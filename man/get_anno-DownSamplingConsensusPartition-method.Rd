\name{get_anno-DownSamplingConsensusPartition-method}
\alias{get_anno,DownSamplingConsensusPartition-method}
\title{
Get annotations
}
\description{
Get annotations
}
\usage{
\S4method{get_anno}{DownSamplingConsensusPartition}(object)
}
\arguments{

  \item{object}{A \code{\link{DownSamplingConsensusPartition-class}} object.}

}
\value{
A data frame if \code{anno} was specified in \code{\link{consensus_partition_by_down_sampling}}, or else \code{NULL}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_ds)
get_anno(golub_cola_ds)
}
