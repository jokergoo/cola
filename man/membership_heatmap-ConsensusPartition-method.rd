\name{membership_heatmap-ConsensusPartition-method}
\alias{membership_heatmap,ConsensusPartition-method}
\alias{membership_heatmap}
\title{
Heatmap of membership of samples in each partition
}
\description{
Heatmap of membership of samples in each partition
}
\usage{
\S4method{membership_heatmap}{ConsensusPartition}(object, k, internal = FALSE,
    anno = get_anno(object), anno_col = get_anno_col(object),
    show_column_names = !internal, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{k}{number of partitions.}
  \item{internal}{used internally.}
  \item{anno}{a data frame with column annotations.}
  \item{anno_col}{colors for the annotations.}
  \item{show_column_names}{whether show column names in the heatmap (which is the column name in the original matrix).}
  \item{...}{other arguments}

}
\details{
Each row in the heatmap is the membership of samples in one partition.

Heatmap is split on rows by \code{top_n}..
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
membership_heatmap(cola_rl["sd", "hclust"], k = 3)
}
