\name{select_partition_number-ConsensusPartition-method}
\alias{select_partition_number,ConsensusPartition-method}
\alias{select_partition_number}
\title{
Several plots for determining the optimized number of partitions
}
\description{
Several plots for determining the optimized number of partitions
}
\usage{
\S4method{select_partition_number}{ConsensusPartition}(object)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}

}
\details{
There are following plots made:

\itemize{
  \item ECDF of the consensus matrix under each k, made by \code{\link{plot_ecdf,ConsensusPartition-method}},
  \item the cophenetic correlation coefficient,
  \item \code{\link{PAC}} score,
  \item mean sihouette score,
  \item the \code{\link{concordance}} for each partition to the consensus partition,
  \item area increase of the area under the ECDF of consensus matrix with increasing k,
  \item Rand index for current k compared to k - 1,
  \item Jaccard coefficient for current k compared to k - 1,
}
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
select_partition_number(cola_rl["sd", "hclust"])
}
