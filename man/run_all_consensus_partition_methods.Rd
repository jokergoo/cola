\name{run_all_consensus_partition_methods}
\alias{run_all_consensus_partition_methods}
\title{
Consensus partitioning for all combinations of methods
}
\description{
Consensus partitioning for all combinations of methods
}
\usage{
run_all_consensus_partition_methods(data,
    top_value_method = all_top_value_methods(),
    partition_method = all_partition_methods(),
    max_k = 6,
    top_n = seq(min(1000, round(nrow(data)*0.1)),
    min(3000, round(nrow(data)*0.3)),
    length.out = 3),
    mc.cores = 1, cores = mc.cores, anno = NULL, anno_col = NULL,
    sample_by = "row", p_sampling = 0.8, partition_repeat = 50,
    scale_rows = NULL, verbose = TRUE, help = cola_opt$help)
}
\arguments{

  \item{data}{A numeric matrix where subgroups are found by columns.}
  \item{top_value_method}{Method which are used to extract top n rows. Allowed methods are in \code{\link{all_top_value_methods}} and can be self-added by \code{\link{register_top_value_methods}}.}
  \item{partition_method}{Method which are used to partition samples.  Allowed methods are in \code{\link{all_partition_methods}} and can be self-added  by \code{\link{register_partition_methods}}.}
  \item{max_k}{Maximal number of subgroups to try. The function will try \code{2:max_k} subgroups.}
  \item{top_n}{Number of rows with top values. The value can be a vector with length > 1. When n > 5000,  the function only randomly sample 5000 rows from top n rows. If \code{top_n} is a vector, paritition will be applied to every values in \code{top_n} and consensus partition is summarized from all partitions.}
  \item{mc.cores}{Number of cores to use. This argument will be removed in future versions.}
  \item{cores}{Number of cores, or a \code{cluster} object returned by \code{\link[parallel]{makeCluster}}.}
  \item{anno}{A data frame with known annotation of columns.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{sample_by}{Should randomly sample the matrix by rows or by columns?}
  \item{p_sampling}{Proportion of the top n rows to sample.}
  \item{partition_repeat}{Number of repeats for the random sampling.}
  \item{scale_rows}{Whether to scale rows. If it is \code{TRUE}, scaling method defined in \code{\link{register_partition_methods}} is used.}
  \item{verbose}{Whether to print messages.}
  \item{help}{Whether to print help messages.}

}
\details{
The function performs consensus partitioning by \code{\link{consensus_partition}} for all combinations of top-value methods and partitioning methods.

It also adjsuts the subgroup labels for all methods and for all k to make them as consistent as possible.
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
}
