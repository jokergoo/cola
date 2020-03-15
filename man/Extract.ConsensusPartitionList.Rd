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

  \item{x}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{i}{Index for top-value methods, character or nummeric.}
  \item{j}{Index for partitioning methods, character or nummeric.}
  \item{drop}{Whether drop class}

}
\details{
For a specific combination of top-value method and partitioning method, you can also
subset by e.g. \code{x['SD:hclust']}.
}
\value{
A \code{\link{ConsensusPartitionList-class}} object or a \code{\link{ConsensusPartition-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
cola_rl[c("SD", "MAD"), c("hclust", "kmeans")]
cola_rl["SD", "kmeans"] # a ConsensusPartition object
cola_rl["SD:kmeans"] # a ConsensusPartition object
cola_rl[["SD:kmeans"]] # a ConsensusPartition object
cola_rl["SD", "kmeans", drop = FALSE] # still a ConsensusPartitionList object
cola_rl["SD:kmeans", drop = FALSE] # still a ConsensusPartitionList object
cola_rl["SD", ]
cola_rl[, "hclust"]
cola_rl[1:2, 1:2]
}
