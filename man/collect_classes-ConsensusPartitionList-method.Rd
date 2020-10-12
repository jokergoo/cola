\name{collect_classes-ConsensusPartitionList-method}
\alias{collect_classes,ConsensusPartitionList-method}
\title{
Collect classes from ConsensusPartitionList object
}
\description{
Collect classes from ConsensusPartitionList object
}
\usage{
\S4method{collect_classes}{ConsensusPartitionList}(object, k, show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    anno = get_anno(object), anno_col = get_anno_col(object),
    simplify = FALSE, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object returned by \code{\link{run_all_consensus_partition_methods}}.}
  \item{k}{Number of subgroups.}
  \item{show_column_names}{Whether to show column names in the heatmap (which is the column name in the original matrix).}
  \item{column_names_gp}{Graphics parameters for column names.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{run_all_consensus_partition_methods}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{simplify}{Internally used.}
  \item{...}{Pass to \code{\link[ComplexHeatmap:draw-HeatmapList-method]{draw,HeatmapList-method}}.}

}
\details{
There are following panels in the plot:

\itemize{
  \item a heatmap showing partitions predicted from all methods where the top annotation is the consensus partition summarized from partitions from all methods, weighted by mean silhouette scores in every single method.
  \item a row barplot annotation showing the mean silhouette scores for different methods.
}

The row clustering is applied on the dissimilarity matrix calculated by \code{\link[clue]{cl_dissimilarity}} with the \code{comembership} method.

The brightness of the color corresponds to the silhouette scores for the consensus partition in each method.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(golub_cola)
collect_classes(golub_cola, k = 3)
}
