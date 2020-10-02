\name{test_to_known_factors-DownSamplingConsensusPartition-method}
\alias{test_to_known_factors,DownSamplingConsensusPartition-method}
\title{
Test correspondance between predicted subgroups and known factors
}
\description{
Test correspondance between predicted subgroups and known factors
}
\usage{
\S4method{test_to_known_factors}{DownSamplingConsensusPartition}(object, k, known = get_anno(object),
    p_cutoff = 0.05, verbose = FALSE)
}
\arguments{

  \item{object}{A \code{\link{DownSamplingConsensusPartition-class}} object.}
  \item{k}{Number of subgroups. It uses all \code{k} if it is not specified.}
  \item{known}{A vector or a data frame with known factors. By default it is the annotation table set in \code{\link{consensus_partition_by_down_sampling}}.}
  \item{p_cutoff}{Cutoff for p-values for the class prediction. Samples with p-value higher than it are omit.}
  \item{verbose}{Whether to print messages.}

}
\details{
The test is performed by \code{\link{test_between_factors}} between the predicted classes and user's annotation table.
}
\value{
A data frame with the following columns:

\itemize{
  \item number of samples used to test after filtered by \code{p_cutoff},
  \item p-values from the tests,
  \item number of subgroups.
}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola_ds)
test_to_known_factors(golub_cola_ds, k = 3)
test_to_known_factors(golub_cola_ds)
}
