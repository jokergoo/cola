\name{top_rows_heatmap-ConsensusPartitionList-method}
\alias{top_rows_heatmap,ConsensusPartitionList-method}
\title{
Heatmap of top rows from different top methods
}
\description{
Heatmap of top rows from different top methods
}
\usage{
\S4method{top_rows_heatmap}{ConsensusPartitionList}(object, top_n = object@list[[1]]@top_n[1], ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{top_n}{number of top rows}
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
