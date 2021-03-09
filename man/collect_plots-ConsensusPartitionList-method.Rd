\name{collect_plots-ConsensusPartitionList-method}
\alias{collect_plots,ConsensusPartitionList-method}
\title{
Collect plots from ConsensusPartitionList object
}
\description{
Collect plots from ConsensusPartitionList object
}
\usage{
\S4method{collect_plots}{ConsensusPartitionList}(object, k = 2, fun = consensus_heatmap,
    top_value_method = object@top_value_method,
    partition_method = object@partition_method,
    verbose = TRUE, mc.cores = 1, cores = mc.cores, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object from \code{\link{run_all_consensus_partition_methods}}.}
  \item{k}{Number of subgroups.}
  \item{fun}{Function used to generate plots. Valid functions are \code{\link{consensus_heatmap}}, \code{\link{plot_ecdf}}, \code{\link{membership_heatmap}}, \code{\link{get_signatures}} and \code{\link{dimension_reduction}}.}
  \item{top_value_method}{A vector of top-value methods.}
  \item{partition_method}{A vector of partitioning methods.}
  \item{verbose}{Whether to print message.}
  \item{mc.cores}{Number of cores. This argument will be removed in figure versions.}
  \item{cores}{Number of cores, or a \code{cluster} object returned by \code{\link[parallel]{makeCluster}}.}
  \item{...}{other Arguments passed to corresponding \code{fun}.}

}
\details{
Plots for all combinations of top-value methods and parittioning methods are arranged in one single page.

This function makes it easy to directly compare results from multiple methods.
}
\value{
No value is returned.
}
\seealso{
\code{\link{collect_plots,ConsensusPartition-method}} collects plots for a single \code{\link{ConsensusPartition-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\donttest{
data(golub_cola)
collect_plots(cola_rl, k = 3)
collect_plots(cola_rl, k = 3, fun = membership_heatmap)
collect_plots(cola_rl, k = 3, fun = get_signatures)
}
}
