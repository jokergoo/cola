\name{membership_heatmap-ConsensusPartition-method}
\alias{membership_heatmap,ConsensusPartition-method}
\alias{membership_heatmap}
\title{
Heatmap of membership in each partition
}
\description{
Heatmap of membership in each partition
}
\usage{
\S4method{membership_heatmap}{ConsensusPartition}(object, k, internal = FALSE,
    anno = get_anno(object), anno_col = get_anno_col(object),
    show_column_names = FALSE, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{k}{Number of subgroups.}
  \item{internal}{Used internally.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{consensus_partition}} or \code{\link{run_all_consensus_partition_methods}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{show_column_names}{Whether show column names in the heatmap (which is the column name in the original matrix).}
  \item{...}{Other arguments}

}
\details{
Each row in the heatmap is the membership in one single partition.

Heatmap is split on rows by \code{top_n}.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
membership_heatmap(cola_rl["SD", "hclust"], k = 3)
}
