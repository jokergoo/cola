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
    top_method = "MAD",
    top_n = seq(min(1000, round(nrow(data)*0.1)),
    min(c(5000, round(nrow(data)*0.5))),
    length.out = 5),
    partition_method = "kmeans",
    k = 2:6, p_sampling = 0.8,
    partition_repeat = 50,
    partition_param = list(),
    known_anno = NULL,
    known_col = NULL,
    scale_rows = NULL,
    verbose = TRUE,
    .env = NULL)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by samples.}
  \item{top_method}{a single top method. Available methods are in \code{\link{all_top_value_methods}}.}
  \item{top_n}{number of rows with top values. When n > 5000, the function only random samples 5000 rows from top n rows.}
  \item{partition_method}{a single partition method. Available ialble methods are in \code{\link{all_partition_methods}}.}
  \item{k}{a list number of partitions.}
  \item{p_sampling}{proportion of the top n rows to sample.}
  \item{partition_repeat}{number of repeats for the random sampling.}
  \item{partition_param}{parameters for the partition method.}
  \item{known_anno}{a data frame with known annotation of samples.}
  \item{known_col}{a list of colors for the annotations in \code{known_anno}.}
  \item{scale_rows}{whether to scale rows. If it is \code{TRUE}, scaling method defined in \code{\link{register_partition_fun}} is used.}
  \item{verbose}{whether print messages}
  \item{.env}{an environment, internally used.}

}
\details{
The function performs analysis by following procedures:

\itemize{
  \item calculate scores for rows by top method and take top n rows
  \item randomly sample \code{p_sampling} rows and perform partitions for \code{partition_repeats} times
  \item collect partitions from all resamplings and calculate consensus partitions
}
}
\value{
A \code{\link{ConsensusPartition-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
