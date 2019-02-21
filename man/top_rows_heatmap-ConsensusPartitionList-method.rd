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
    scale_rows = object@list[[1]]@scale_rows, ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{top_n}{number of top rows.}
  \item{scale_rows}{wether scale rows. }
  \item{...}{pass to \code{\link{top_rows_heatmap,matrix-method}}}

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
# There is no example
NULL
}
