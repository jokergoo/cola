\name{top_rows_heatmap-ConsensusPartition-method}
\alias{top_rows_heatmap,ConsensusPartition-method}
\title{
Heatmap of top rows
}
\description{
Heatmap of top rows
}
\usage{
\S4method{top_rows_heatmap}{ConsensusPartition}(object, top_n = min(object@top_n), k = NULL,
    anno = get_anno(object), anno_col = get_anno_col(object),
    scale_rows = object@scale_rows, ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{top_n}{Number of top rows.}
  \item{k}{-k}
  \item{anno}{A data frame of annotations.}
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
top_rows_heatmap(golub_cola["ATC:skmeans"])
}
}
