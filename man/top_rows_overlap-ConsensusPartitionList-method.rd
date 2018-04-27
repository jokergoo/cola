\name{top_rows_overlap-ConsensusPartitionList-method}
\alias{top_rows_overlap,ConsensusPartitionList-method}
\title{
Overlap of top rows from different top value methods
}
\description{
Overlap of top rows from different top value methods
}
\usage{
\S4method{top_rows_overlap}{ConsensusPartitionList}(object, top_n = min(object@list[[1]]@top_n),
    method = c("venn", "venneuler", "correspondance"), ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{top_n}{number of top rows.}
  \item{method}{\code{venn}: use Venn diagram; \code{venneuler}: use Venn Euler diagram; \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{...}{additional arguments passed to \code{\link{venn_euler}} or \code{\link{correspond_between_rankings}}.}

}
\value{
No value is returned.
}
\seealso{
\code{\link{top_rows_overlap,list-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
data(cola_rl)
top_rows_overlap(cola_rl, method = "venn")
top_rows_overlap(cola_rl, method = "correspondance")
}
