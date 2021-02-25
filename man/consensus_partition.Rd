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
    top_value_method = "ATC",
    top_n = seq(min(1000, round(nrow(data)*0.1)),
    min(3000, round(nrow(data)*0.3)),
    length.out = 3),
    partition_method = "skmeans",
    max_k = 6,
    sample_by = "row",
    p_sampling = 0.8,
    partition_repeat = 50,
    partition_param = list(),
    anno = NULL,
    anno_col = NULL,
    scale_rows = NULL,
    verbose = TRUE,
    mc.cores = 1, cores = mc.cores,
    prefix = "",
    .env = NULL,
    help = cola_opt$help)
}
\arguments{

  \item{data}{A numeric matrix where subgroups are found by columns.}
  \item{top_value_method}{A single top-value method. Available methods are in \code{\link{all_top_value_methods}}. Use \code{\link{register_top_value_methods}} to add a new top-value method.}
  \item{top_n}{Number of rows with top values. The value can be a vector with length > 1. When n > 5000,  the function only randomly sample 5000 rows from top n rows. If \code{top_n} is a vector, paritition will be applied to every values in \code{top_n} and consensus partition is summarized from all partitions.}
  \item{partition_method}{A single partitioning method. Available methods are in \code{\link{all_partition_methods}}. Use \code{\link{register_partition_methods}} to add a new partition method.}
  \item{max_k}{Maximal number of subgroups to try. The function will try for \code{2:max_k} subgroups}
  \item{sample_by}{Should randomly sample the matrix by rows or by columns?}
  \item{p_sampling}{Proportion of the submatrix which contains the top n rows to sample.}
  \item{partition_repeat}{Number of repeats for the random sampling.}
  \item{partition_param}{Parameters for the partition method which are passed to \code{...} in a registered partitioning method. See \code{\link{register_partition_methods}} for detail.}
  \item{anno}{A data frame with known annotation of samples. The annotations will be plotted in heatmaps and the correlation to predicted subgroups will be tested.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{scale_rows}{Whether to scale rows. If it is \code{TRUE}, scaling method defined in \code{\link{register_partition_methods}} is used.}
  \item{verbose}{Whether print messages.}
  \item{mc.cores}{Multiple cores to use. This argument will be removed in future versions.}
  \item{cores}{Multiple cores to use.}
  \item{prefix}{Internally used.}
  \item{.env}{An environment, internally used.}
  \item{help}{Whether to print help messages.}

}
\details{
The function performs analysis in following steps:

\itemize{
  \item calculate scores for rows by top-value method,
  \item for each top_n value, take top n rows,
  \item randomly sample \code{p_sampling} rows from the top_n-row matrix and perform partitioning for \code{partition_repeats} times,
  \item collect partitions from all individual partitions and summarize a consensus partition.
}
}
\value{
A \code{\link{ConsensusPartition-class}} object. Simply type object in the interactive R session
to see which functions can be applied on it.
}
\seealso{
\code{\link{run_all_consensus_partition_methods}} runs consensus partitioning with multiple top-value methods
and multiple partitioning methods.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
set.seed(123)
m = cbind(rbind(matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
                matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
                matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
          rbind(matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20),
                matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20),
                matrix(rnorm(20*20, mean = 0,   sd = 0.5), nr = 20)),
          rbind(matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
                matrix(rnorm(20*20, mean = 0.5, sd = 0.5), nr = 20),
                matrix(rnorm(20*20, mean = 1,   sd = 0.5), nr = 20))
         ) + matrix(rnorm(60*60, sd = 0.5), nr = 60)
res = consensus_partition(m, partition_repeat = 10, top_n = c(10, 20, 50))
res
}
