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
    top_method = all_top_value_methods()[1],
    top_n = seq(min(1000, round(nrow(data)*0.1)),
    min(c(5000, round(nrow(data)*0.5))),
    length.out = 5),
    partition_method = all_partition_methods()[1],
    k = 2:6, p_sampling = 0.8,
    partition_repeat = 50,
    partition_param = list(),
    known_anno = NULL,
    known_col = NULL,
    scale_rows = NULL,
    column_index = seq_len(ncol(data)),
    .env = NULL)
}
\arguments{

  \item{data}{a numeric matrix where subgroups are found by columns.}
  \item{top_method}{a single top method. Available methods are in \code{\link{all_top_value_methods}}.}
  \item{top_n}{number of rows with top values. When n > 5000, the function only random samples 5000 rows from top n rows.}
  \item{partition_method}{a single partition method. Available ialble methods are in \code{\link{all_partition_methods}}.}
  \item{k}{a list number of partitions.}
  \item{p_sampling}{proportion of the top n rows to sample.}
  \item{partition_repeat}{number of repeats for the random sampling.}
  \item{partition_param}{parameters for the partition method.}
  \item{known_anno}{a data frame with known annotation of columns.}
  \item{known_col}{a list of colors for the annotations in \code{known_anno}.}
  \item{scale_rows}{whether to scale rows.}
  \item{column_index}{a subst of columns in the matrix, only for internal use.}
  \item{.env}{an environment, internally used.}

}
\details{
The function performs analysis by following procedures:

\itemize{
  \item calculate scores for rows by top method and take top n rows
  \item randomly sample \code{p_sampling} rows and perform partitions for \code{partition_repeats} times
  \item collect partitions from all randomizations can calculate consensus clusters
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
