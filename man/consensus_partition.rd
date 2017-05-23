\name{consensus_partition}
\alias{consensus_partition}
\title{
Consensus partition
}
\description{
Consensus partition
}
\usage{
consensus_partition(data,
    top_method = ALL_TOP_VALUE_METHOD()[1],
    top_n = seq(min(2000, round(nrow(data)*0.2)), min(c(6000, round(nrow(data)*0.6))), length.out = 5),
    partition_method = ALL_PARTITION_METHOD()[1],
    k = 2:6, p_sampling = 0.8,
    partition_repeat = 50,
    partition_param = list(),
    known_anno = NULL,
    known_col = NULL,
    .env)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_method}{a single top method. Avaialble methods are in \code{\link{ALL_TOP_VALUE_METHOD}}.}
  \item{top_n}{number of rows with top values.}
  \item{partition_method}{a single partition method. Avaialble methods are in \code{\link{ALL_PARTITION_METHOD}}.}
  \item{k}{number of partitions. The value is a vector.}
  \item{p_sampling}{proportion of the top n rows to sample.}
  \item{partition_repeat}{number of repeats for the random sampling.}
  \item{partition_param}{parameters for the partition method.}
  \item{known_anno}{a data frame with known annotation of columns.}
  \item{known_col}{a list of colors for the annotations in \code{known_anno}.}
  \item{.env}{an environment, internally used.}

}
\value{
A \code{\link{ConsensusPartition-class}} object.
}
\examples{
# There is no example
NULL

}
