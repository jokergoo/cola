\name{[[.ConsensusPartitionList}
\alias{[[.ConsensusPartitionList}
\alias{ExtractExtract.ConsensusPartitionList}
\title{
Subset a ConsensusPartitionList object
}
\description{
Subset a ConsensusPartitionList object
}
\usage{
\method{[[}{ConsensusPartitionList}(x, i)
}
\arguments{

  \item{x}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{i}{character index for combination of top-value methods and partition method.}

}
\value{
A \code{\link{ConsensusPartition-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
cola_rl[["sd:MAD"]]
}
