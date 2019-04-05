\name{collect_stats-ConsensusPartitionList-method}
\alias{collect_stats,ConsensusPartitionList-method}
\title{
Draw and compare statistics for multiple methods
}
\description{
Draw and compare statistics for multiple methods
}
\usage{
\S4method{collect_stats}{ConsensusPartitionList}(object, k, layout_nrow = 2, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{k}{number of partitions}
  \item{layout_nrow}{number of rows in the layout}
  \item{...}{other arguments}

}
\details{
It draws heatmaps for statistics for multiple methods in parallel, so that users can compare which combination
of methods gives the best results with given the number of partitions.
}
\examples{
data(cola_rl)
collect_stats(cola_rl, k = 3)
}
