\name{[.ConsensusPartitionList}
\alias{[.ConsensusPartitionList}
\alias{Extract.ConsensusPartitionList}
\title{
Subset a ConsensusPartitionList object
}
\description{
Subset a ConsensusPartitionList object
}
\usage{
\method{[}{ConsensusPartitionList}(x, i, j, drop = TRUE)
}
\arguments{

  \item{x}{a \code{\link{ConsensusPartition-class}} object.}
  \item{i}{index for top value methods, character or nummeric.}
  \item{j}{index for partition methods, character or nummeric.}
  \item{drop}{whether drop class}

}
\value{
A \code{\link{ConsensusPartitionList-class}} object or a \code{\link{ConsensusPartition-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
cola_rl[c("sd", "MAD"), c("hclust", "kmeans")]
cola_rl["sd", "kmeans"]
cola_rl["sd", "kmeans", drop = FALSE]
cola_rl["sd", ]
cola_rl[, "hclust"]
cola_rl[1:2, 1:2]
}
