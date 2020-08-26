\name{plot_ecdf-ConsensusPartition-method}
\alias{plot_ecdf,ConsensusPartition-method}
\alias{plot_ecdf}
\title{
Plot the empirical cumulative distribution (eCDF) curve of the consensus matrix
}
\description{
Plot the empirical cumulative distribution (eCDF) curve of the consensus matrix
}
\usage{
\S4method{plot_ecdf}{ConsensusPartition}(object, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{...}{Other arguments.}

}
\details{
It plots eCDF curve for each k.

This function is mainly used in \code{\link{collect_plots}} and \code{\link{select_partition_number}} functions.
}
\value{
No value is returned.
}
\seealso{
See \code{\link[stats]{ecdf}} for a detailed explanation of the empirical cumulative distribution function.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
plot_ecdf(golub_cola["ATC", "skmeans"])
}
