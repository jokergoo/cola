\name{collect_plots-ConsensusPartition-method}
\alias{collect_plots,ConsensusPartition-method}
\title{
Collect plots from ConsensusPartition object
}
\description{
Collect plots from ConsensusPartition object
}
\usage{
\S4method{collect_plots}{ConsensusPartition}(object, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{...}{other arguments}

}
\details{
Plots by \code{\link{plot_ecdf,ConsensusPartition-method}}, \code{\link{collect_classes,ConsensusPartition-method}}, \code{\link{consensus_heatmap,ConsensusPartition-method}}, \code{\link{membership_heatmap,ConsensusPartition-method}} 
and \code{\link{get_signatures,ConsensusPartition-method}} are arranged in one single page, for all avaialble k.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
