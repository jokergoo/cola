\name{cola_report-ConsensusPartition-method}
\alias{cola_report,ConsensusPartition-method}
\title{
Make report from the ConsensusPartition object
}
\description{
Make report from the ConsensusPartition object
}
\usage{
\S4method{cola_report}{ConsensusPartition}(object, output_dir)
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartition-class}} object.}
  \item{output_dir}{the output directory where the report is put.}

}
\details{
Please generate report on the \code{\link{ConsensusPartitionList-class}} object directly.

If you want to make report only for one single method, you can subset the 
\code{\link{ConsensusPartitionList-class}} object and then call \code{cola_report}, e.g.

  \preformatted{
    cola_report(res_list["sd", "hclust", drop = FALSE], output_dir = ...)  }
}
\value{
No value is returned.
}
\seealso{
\code{\link{cola_report,ConsensusPartitionList-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
