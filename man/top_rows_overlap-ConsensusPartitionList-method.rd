\name{top_rows_overlap-ConsensusPartitionList-method}
\alias{top_rows_overlap,ConsensusPartitionList-method}
\title{
Overlap of top rows from different top methods
}
\description{
Overlap of top rows from different top methods
}
\usage{
\S4method{top_rows_overlap}{ConsensusPartitionList}(object, top_n = round(0.25*length(all_value_list[[1]])),
    type = c("venn", "correspondance"), ...)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{top_n}{number of top rows}
  \item{type}{\code{venn}: use Venn Euler digram; \code{correspondance}: use \code{\link{correspond_between_rankings}}.}
  \item{...}{additional arguments passed to \code{\link{venn_euler}} or \code{\link{correspond_between_rankings}}}

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
# There is no example
NULL

}
