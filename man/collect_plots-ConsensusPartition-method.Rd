\name{collect_plots-ConsensusPartition-method}
\alias{collect_plots,ConsensusPartition-method}
\title{
Collect plots from ConsensusPartition object
}
\description{
Collect plots from ConsensusPartition object
}
\usage{
\S4method{collect_plots}{ConsensusPartition}(object, verbose = TRUE)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{verbose}{Whether print messages.}

}
\details{
Plots by \code{\link{plot_ecdf}}, \code{\link{collect_classes,ConsensusPartition-method}}, \code{\link{consensus_heatmap}}, \code{\link{membership_heatmap}} 
and \code{\link{get_signatures}} are arranged in one single page, for all avaialble k.
}
\value{
No value is returned.
}
\seealso{
\code{\link{collect_plots,ConsensusPartitionList-method}} collects plots for the \code{\link{ConsensusPartitionList-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rl)
collect_plots(cola_rl["SD", "kmeans"])
}
}
