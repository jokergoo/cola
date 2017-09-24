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
    anno_col = if(missing(anno)) object@known_col else NULL,
    show_row_names = FALSE, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}
  \item{show_legend}{whether show heatmap and annotation legends.}
  \item{anno}{a data frame with column annotations}
  \item{anno_col}{colors for the annotations}
  \item{show_row_names}{whether plot row names on the consensus heatmap}
  \item{...}{other arguments}

}
\details{
There are following heatmaps from left to right:

\itemize{
  \item probability of the column to stay in the subgroup
  \item silhouette values which measure the distance for an item to the second closest subgroups
  \item predicted classes
  \item consensus matrix
  \item more annotations if provided as \code{anno}
}
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
