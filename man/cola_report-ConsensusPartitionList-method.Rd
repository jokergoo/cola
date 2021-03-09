\name{cola_report-ConsensusPartitionList-method}
\alias{cola_report,ConsensusPartitionList-method}
\title{
Make HTML report from the ConsensusPartitionList object
}
\description{
Make HTML report from the ConsensusPartitionList object
}
\usage{
\S4method{cola_report}{ConsensusPartitionList}(object, output_dir = getwd(), mc.cores = 1, cores = mc.cores,
    title = "cola Report for Consensus Partitioning", env = parent.frame())
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{output_dir}{The output directory where the report is saved.}
  \item{mc.cores}{Multiple cores to use. This argument will be removed in future versions.}
  \item{cores}{Number of cores, or a \code{cluster} object returned by \code{\link[parallel]{makeCluster}}.}
  \item{title}{Title of the report.}
  \item{env}{Where the objects in the report are found, internally used.}

}
\details{
The \code{\link{ConsensusPartitionList-class}} object contains results for all combinations of top-value methods and partitioning methods.
This function generates a HTML report which contains all plots and tables for every combination
of method.

The report generation may take a while because it generates A LOT of heatmaps.

Examples of reports can be found at \url{https://jokergoo.github.io/cola_collection/} .
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
if(FALSE) {
# the following code is runnable
data(golub_cola)
cola_report(golub_cola[c("SD", "MAD"), c("hclust", "skmeans")], output_dir = "~/test_cola_cl_report")
}
}
