\name{run_all_consensus_partition_methods}
\alias{run_all_consensus_partition_methods}
\title{
Run subgroup classification for all methods
}
\description{
Run subgroup classification for all methods
}
\usage{
run_all_consensus_partition_methods(data, top_method = all_top_value_methods(),
    partition_method = all_partition_methods(),
    mc.cores = 1, ...)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_method}{method which are used to extract top n rows. Allowed methods are in \code{\link{all_top_value_methods}} and can be self-added by \code{\link{register_top_value_fun}}.}
  \item{partition_method}{method which are used to do partition on data columns.  Allowed methods are in \code{\link{all_partition_methods}} and can be self-added  by \code{\link{register_partition_fun}}.}
  \item{mc.cores}{number of cores to use`.}
  \item{...}{other arguments passed to \code{\link{consensus_partition}}.}

}
\details{
The function run consensus partitions for all combination of top methods and parittion methods.
}
\value{
a \code{\link{ConsensusPartitionList-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
