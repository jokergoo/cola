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

  \item{object}{a \code{\link{ConsensusPartitionList-class}} object.}
  \item{output_dir}{the output directory where the report is put.}
  \item{env}{where the objects in the report are found, internally used.}

}
\details{
The \code{\link{ConsensusPartitionList-class}} object contains results for all top value methods and all partition methods.
This function generates a HTML report which contains all plots for every combination
of top value method and partition method. The report generation may take a while.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
\dontrun{
data(cola_rl)
cola_report(cola_rl[c("sd", "MAD"), c("hclust", "skmeans")], output_dir = "~/test")
}
}
