\name{run_all_consensus_partition_methods}
\alias{run_all_consensus_partition_methods}
\title{
Run subgroup classification from all methods
}
\description{
Run subgroup classification from all methods
}
\usage{
run_all_consensus_partition_methods(data, top_method = all_top_value_methods(),
    partition_method = all_partition_methods(), k = 2:6,
    mc.cores = 1, known_anno = NULL, known_col = NULL, ...)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by samples}
  \item{top_method}{method which are used to extract top n rows. Allowed methods are in \code{\link{all_top_value_methods}} and can be self-added by \code{\link{register_top_value_fun}}.}
  \item{partition_method}{method which are used to do partition on samples.  Allowed methods are in \code{\link{all_partition_methods}} and can be self-added  by \code{\link{register_partition_fun}}.}
  \item{k}{a list number of partitions.}
  \item{mc.cores}{number of cores to use.}
  \item{known_anno}{a data frame with known annotation of samples}
  \item{known_col}{a list of colors for the annotations in \code{known_anno}.}
  \item{...}{other arguments passed to \code{\link{consensus_partition}}.}

}
\details{
The function run consensus partitions for all combinations of top methods and parittion methods.
}
\value{
A \code{\link{ConsensusPartitionList-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
