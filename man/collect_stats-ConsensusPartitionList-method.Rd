\name{collect_stats-ConsensusPartitionList-method}
\alias{collect_stats,ConsensusPartitionList-method}
\title{
Draw and compare statistics for multiple methods
}
\description{
Draw and compare statistics for multiple methods
}
\usage{
\S4method{collect_stats}{ConsensusPartitionList}(object, k, layout_nrow = 2, all_stats = FALSE, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{k}{Number of subgroups.}
  \item{layout_nrow}{Number of rows in the layout}
  \item{all_stats}{Whether to show all statistics that were calculated. Used internally.}
  \item{...}{Other arguments}

}
\details{
It draws heatmaps for statistics for multiple methods in parallel, so that users can compare which combination
of methods gives the best results with given the number of subgroups.
}
\examples{
data(cola_rl)
collect_stats(cola_rl, k = 3)
}
