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
data(golub_cola)
golub_cola[c("SD", "MAD"), c("hclust", "kmeans")]
golub_cola["SD", "kmeans"] # a ConsensusPartition object
golub_cola["SD:kmeans"] # a ConsensusPartition object
golub_cola[["SD:kmeans"]] # a ConsensusPartition object
golub_cola["SD", "kmeans", drop = FALSE] # still a ConsensusPartitionList object
golub_cola["SD:kmeans", drop = FALSE] # still a ConsensusPartitionList object
golub_cola["SD", ]
golub_cola[, "hclust"]
golub_cola[1:2, 1:2]
}
