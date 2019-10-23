\name{top_rows_overlap-ConsensusPartitionList-method}
\alias{top_rows_overlap,ConsensusPartitionList-method}
\title{
Overlap of top rows from different top-value methods
}
\description{
Overlap of top rows from different top-value methods
}
\usage{
\S4method{top_rows_overlap}{ConsensusPartitionList}(object, top_n = min(object@list[[1]]@top_n),
    method = c("euler", "venn", "correspondance"), ...)
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartitionList-class}} object.}
  \item{top_n}{Number of top rows.}
  \item{method}{\code{euler}: plot Euler diagram by \code{\link[eulerr]{euler}}; \code{venn}: plot Venn diagram by \code{\link[gplots]{venn}};  \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{...}{Additional arguments passed to \code{\link[eulerr]{plot.euler}} or \code{\link{correspond_between_rankings}}.}

}
\value{
No value is returned.
}
\seealso{
\code{\link{top_elements_overlap}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
top_rows_overlap(cola_rl, method = "venn")
top_rows_overlap(cola_rl, method = "correspondance")
}
