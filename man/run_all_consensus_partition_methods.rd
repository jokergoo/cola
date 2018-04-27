\name{run_all_consensus_partition_methods}
\alias{run_all_consensus_partition_methods}
\title{
Run subgroup classification from all methods
}
\description{
Run subgroup classification from all methods
}
\usage{
run_all_consensus_partition_methods(data,
    top_value_method = all_top_value_methods(),
    partition_method = all_partition_methods(), max_k = 6,
    mc.cores = 1, anno = NULL, anno_col = NULL, ...)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_value_method}{method which are used to extract top n rows. Allowed methods are in \code{\link{all_top_value_methods}} and can be self-added by \code{\link{register_top_value_method}}.}
  \item{partition_method}{method which are used to do partition on samples.  Allowed methods are in \code{\link{all_partition_methods}} and can be self-added  by \code{\link{register_partition_method}}.}
  \item{max_k}{maximum number of partitions to try. It will use \code{2:max_k} partitions.}
  \item{mc.cores}{number of cores to use.}
  \item{anno}{a data frame with known annotation of samples}
  \item{anno_col}{a list of colors for the annotations in \code{known_anno}.}
  \item{...}{other arguments passed to \code{\link{consensus_partition}}.}

}
\details{
The function run consensus partitions for all combinations of top methods and parittion methods.
}
\value{
A \code{\link{ConsensusPartitionList-class}} object. Simply type the name of the object in the R interactive session
to see which functions can be applied on it.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
set.seed(123)
m = cbind(rbind(matrix(rnorm(20*20, mean = 1), nr = 20),
			    matrix(rnorm(20*20, mean = -1), nr = 20)),
	rbind(matrix(rnorm(20*20, mean = -1), nr = 20),
			    matrix(rnorm(20*20, mean = 1), nr = 20))
) + matrix(rnorm(40*40), nr = 40)
rl = run_all_consensus_partition_methods(data = m, top_n = c(20, 30, 40))
}
data(cola_rl)
cola_rl
}
