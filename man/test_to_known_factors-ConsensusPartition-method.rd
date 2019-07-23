\name{test_to_known_factors-ConsensusPartition-method}
\alias{test_to_known_factors,ConsensusPartition-method}
\title{
Test correspondance between predicted classes and known factors
}
\description{
Test correspondance between predicted classes and known factors
}
\usage{
\S4method{test_to_known_factors}{ConsensusPartition}(object, k, known = get_anno(object),
    silhouette_cutoff = 0.5, verbose = FALSE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of partitions. It uses all \code{k} if it is not set.}
  \item{known}{A vector or a data frame with known factors. By default it is the annotation table set in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.}
  \item{silhouette_cutoff}{Cutoff for sihouette scores. Samples with value less than it are omit.}
  \item{verbose}{Whether to print messages.}

}
\value{
A data frame with columns:

\itemize{
  \item number of samples used to test after filtered by \code{silhouette_cutoff}
  \item p-values from the tests
  \item number of partitions
}
}
\seealso{
\code{\link{test_between_factors}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
test_to_known_factors(cola_rl[1, 1], known = 1:40)
}
