\name{run_all_consensus_partition_methods}
\alias{run_all_consensus_partition_methods}
\title{
Consensus partition for all combinations of methods
}
\description{
Consensus partition for all combinations of methods
}
\usage{
run_all_consensus_partition_methods(data,
    top_value_method = all_top_value_methods(),
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
  \item{max_k}{maximal number of partitions to try. The function will try \code{2:max_k} partitions.}
  \item{top_n}{number of rows with top values. The value can be a vector with length > 1. When n > 5000,  the function only randomly sample 5000 rows from top n rows. If \code{top_n} is a vector, paritition will be applied to every values in \code{top_n} and consensus partition is summarized from all partitions.}
  \item{mc.cores}{number of cores to use.}
  \item{anno}{a data frame with known annotation of columns.}
  \item{anno_col}{a list of colors for the annotations in \code{anno}.}
  \item{p_sampling}{proportion of the top n rows to sample.}
  \item{partition_repeat}{number of repeats for the random sampling.}
  \item{scale_rows}{whether to scale rows. If it is \code{TRUE}, scaling method defined in \code{\link{register_partition_method}} is used.}

}
\details{
The function runs consensus partitions by \code{\link{consensus_partition}} for all combinations of top methods and parittion methods.

It also adjsuts the class IDs for all methods and for all k to make them as consistent as possible.
}
\value{
A \code{\link{ConsensusPartitionList-class}} object. Simply type object in the interactive R session
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
