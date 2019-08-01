\name{collect_plots-HierarchicalPartition-method}
\alias{collect_plots,HierarchicalPartition-method}
\title{
Collect plots from HierarchicalPartition object
}
\description{
Collect plots from HierarchicalPartition object
}
\usage{
\S4method{collect_plots}{HierarchicalPartition}(object, depth = max_depth(object),
    fun = consensus_heatmap, verbose = TRUE, mc.cores = 1, ...)
}
\arguments{

  \item{object}{A \code{\link{HierarchicalPartition-class}} object.}
  \item{depth}{Depth in the hierarchy.}
  \item{fun}{Function used to generate plots. Valid functions are \code{\link{consensus_heatmap}}, \code{\link{plot_ecdf}}, \code{\link{membership_heatmap}}, \code{\link{get_signatures}} and \code{\link{dimension_reduction}}.}
  \item{verbose}{Whether to print message.}
  \item{mc.cores}{Number of cores. On OSX it is enforced to be 1.}
  \item{...}{other Arguments passed to corresponding \code{fun}.}

}
\details{
The hierarchy represents as a circular dendrogram where plots are on the nodes.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rh)
collect_plots(cola_rh)
\dontrun{
collect_plots(cola_rh, fun = membership_heatmap)
collect_plots(cola_rh, fun = get_signatures)
}
}
