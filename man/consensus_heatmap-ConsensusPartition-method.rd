\name{consensus_heatmap-ConsensusPartition-method}
\alias{consensus_heatmap,ConsensusPartition-method}
\alias{consensus_heatmap}
\title{
Heatmap for the consensus matrix
}
\description{
Heatmap for the consensus matrix
}
\usage{
\S4method{consensus_heatmap}{ConsensusPartition}(object, k, show_legend = TRUE,
    anno = object@known_anno,
    anno_col = if(missing(anno)) object@known_col else NULL, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}
  \item{show_legend}{whether show heatmap and annotation legends.}
  \item{anno}{a data frame with column annotations}
  \item{anno_col}{colors for the annotations}
  \item{...}{other arguments}

}
\examples{
# There is no example
NULL

}
