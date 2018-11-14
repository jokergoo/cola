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
    top_value_method = setdiff(all_top_value_methods(), "som"),
    partition_method = all_partition_methods(),
    max_k = 6,
    top_n = seq(min(1000, round(nrow(data)*0.1)),
    min(5000, round(nrow(data)*0.5)),
    length.out = 5),
    mc.cores = 1, anno = NULL, anno_col = NULL,
    p_sampling = 0.8, partition_repeat = 50, scale_rows = NULL)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_value_method}{method which are used to extract top n rows. Allowed methods are in \code{\link{all_top_value_methods}} and can be self-added by \code{\link{register_top_value_method}}.}
  \item{partition_method}{method which are used to do partition on samples.  Allowed methods are in \code{\link{all_partition_methods}} and can be self-added  by \code{\link{register_partition_method}}.}
  \item{max_k}{maximum number of partitions to try. It will use \code{2:max_k} partitions.}
  \item{top_n}{number of rows with top values. The value can be a vector with length > 1. When n > 5000, the function only randomly sample 5000 rows from top n rows.}
  \item{mc.cores}{number of cores to use.}
  \item{anno}{a data frame with known annotation of samples}
  \item{anno_col}{a list of colors for the annotations in \code{known_anno}.}
  \item{p_sampling}{proportion of the top n rows to sample.}
  \item{partition_repeat}{number of repeats for the random sampling.}
  \item{scale_rows}{whether to scale rows. If it is \code{TRUE}, scaling method defined in \code{\link{register_partition_method}} is used.}

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
