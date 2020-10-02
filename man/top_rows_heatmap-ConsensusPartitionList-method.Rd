\name{top_rows_heatmap-ConsensusPartitionList-method}
\alias{top_rows_heatmap,ConsensusPartitionList-method}
\title{
Heatmap of top rows from different top-value methods
}
\description{
Heatmap of top rows from different top-value methods
}
\usage{
\S4method{top_rows_heatmap}{ConsensusPartitionList}(object, top_n = min(object@list[[1]]@top_n),
    anno = get_anno(object), anno_col = get_anno_col(object),
    scale_rows = object@list[[1]]@scale_rows, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{top_n}{Number of top rows.}
  \item{anno}{A data frame of annotations for the original matrix columns.  By default it uses the annotations specified in \code{\link{run_all_consensus_partition_methods}}.}
  \item{anno_col}{A list of colors (color is defined as a named vector) for the annotations. If \code{anno} is a data frame, \code{anno_col} should be a named list where names correspond to the column names in \code{anno}.}
  \item{scale_rows}{Wether to scale rows. }
  \item{...}{Pass to \code{\link{top_rows_heatmap,matrix-method}}}

}
\value{
No value is returned.
}
\seealso{
\code{\link{top_rows_heatmap,matrix-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\donttest{
data(golub_cola)
top_rows_heatmap(golub_cola)
}
}
