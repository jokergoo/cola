\name{run_all_consensus_partition_methods}
\alias{run_all_consensus_partition_methods}
\title{
Run subgroup classification in a batch
}
\description{
Run subgroup classification in a batch
}
\usage{
run_all_consensus_partition_methods(data, top_method = ALL_TOP_VALUE_METHOD(),
    partition_method = ALL_PARTITION_METHOD(),
    mc.cores = 1, get_signatures = FALSE, ...)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_method}{method which are used to extract top n rows. Allowed methods are in \code{\link{ALL_TOP_VALUE_METHOD}} and can be self-added by \code{\link{register_top_value_fun}}.}
  \item{partition_method}{method which are used to do partition on data columns.  Allowed methods are in \code{\link{ALL_PARTITION_METHOD}} and can be self-added  by \code{\link{register_partition_fun}}.}
  \item{mc.cores}{number of cores to use`.}
  \item{get_signatures}{whether to run \code{\link{get_signatures}} for each partition.}
  \item{...}{other arguments passed to \code{\link{consensus_partition}}.}

}
\value{
a \code{\link{ConsensusPartitionList-class}} object.
}
\examples{
# There is no example
NULL

}
