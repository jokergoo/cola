\name{cola_report-ConsensusPartitionList-method}
\alias{cola_report,ConsensusPartitionList-method}
\title{
Make report for the ConsensusPartitionList object
}
\description{
Make report for the ConsensusPartitionList object
}
\usage{
\S4method{cola_report}{ConsensusPartitionList}(object, output_dir = getwd(), env = parent.frame())
}
\arguments{

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object}
  \item{output_dir}{the output directory where put the report}
  \item{env}{where the objects in the report are found, internally used}

}
\details{
The \code{\link{ConsensusPartitionListclass}} object contains results for all top methods and all partition methods.
This function generates a HTML report which contains all plots for every combination
of top method and partition method.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
