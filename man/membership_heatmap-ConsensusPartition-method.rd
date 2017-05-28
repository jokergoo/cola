\name{membership_heatmap-ConsensusPartition-method}
\alias{membership_heatmap,ConsensusPartition-method}
\alias{membership_heatmap}
\title{
Heatmap of membership of columns in each random sampling
}
\description{
Heatmap of membership of columns in each random sampling
}
\usage{
\S4method{membership_heatmap}{ConsensusPartition}(object, k, show_legend = TRUE,
    anno = object@known_anno,
    anno_col = if(missing(anno)) object@known_col else NULL,
    show_column_names = TRUE, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object}
  \item{k}{number of partitions.}
  \item{show_legend}{whether show heatmap and annotation legends.}
  \item{anno}{a data frame with column annotations}
  \item{anno_col}{colors for the annotations}
  \item{show_column_names}{whether show column names in the heatmap}
  \item{...}{other arguments}

}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
`
}
\examples{
# There is no example
NULL

}
