\name{suggest_best_k-ConsensusPartitionList-method}
\alias{suggest_best_k,ConsensusPartitionList-method}
\title{
Suggest the best number of subgroups
}
\description{
Suggest the best number of subgroups
}
\usage{
\S4method{suggest_best_k}{ConsensusPartitionList}(object, jaccard_index_cutoff = 0.95)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{jaccard_index_cutoff}{The cutoff for Jaccard index for comparing to previous k.}

}
\details{
It basically gives the best k for each combination of top-value method and partitioning method by calling \code{\link{suggest_best_k,ConsensusPartition-method}}.

1-PAC score higher than 0.95 is treated as very stable partition (marked by \code{**}) and higher than 0.9 is treated as stable partition (marked by \code{*}).
}
\value{
A data frame with the best k and other statistics for each combination of methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
suggest_best_k(golub_cola)
}
