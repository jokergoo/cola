\name{cola_report-ConsensusPartition-method}
\alias{cola_report,ConsensusPartition-method}
\title{
Make HTML report from the ConsensusPartition object
}
\description{
Make HTML report from the ConsensusPartition object
}
\usage{
\S4method{cola_report}{ConsensusPartition}(object, output_dir = getwd(),
    title = qq("cola Report for Consensus Partitioning (@{object@top_value_method}:@{object@partition_method})"),
    env = parent.frame())
}
\arguments{

  \item{object}{A \code{\link{ConsensusPartition-class}} object.}
  \item{output_dir}{The output directory where the report is saved.}
  \item{title}{Title of the report.}
  \item{env}{Where the objects in the report are found, internally used.}

}
\details{
It generates report for a specific combination of top-value method and partition method.
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
